{-# LANGUAGE TypeFamilies, FlexibleContexts, ScopedTypeVariables #-}

module Main where

import Control.Applicative
import Control.Monad.Reader
import qualified Data.ByteString.Char8 as BS
import Data.List
import Data.Maybe
import Numeric
import System.Console.GetOpt
import System.Environment
import System.IO

import qualified Data.Attoparsec.Char8 as AP
import qualified Data.Vector.Unboxed as U

import Statistics.Quantile

import qualified Bio.SamTools.BamIndex as BamIndex
import qualified Bio.SamTools.FaIdx as FaIdx
import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.Position as Pos
import Bio.SeqLoc.Strand
import Bio.SeqLoc.Transcript

import Bio.RiboSeq.BamFile
import Bio.RiboSeq.CodonAssignment
import Bio.RiboSeq.Translation

data Conf = Conf { confBed :: !FilePath
                 , confOut :: !FilePath
                 , confASite :: !FilePath
                 , confFasta :: !(Maybe FilePath)
                 , confMinMedian :: !Double
                 , confMinPause :: !Double
                 }
          deriving (Show)
        
defaultMinMedian :: Double                   
defaultMinMedian = 3.0 - 1e-6

defaultMinPause :: Double
defaultMinPause = 30.0 - 1e-6
                   
main :: IO ()
main = getArgs >>= handleOpt . getOpt RequireOrder optDescrs
    where handleOpt (_,    _,         errs@(_:_)) = usage (unlines errs)
          handleOpt (args, [bam], []) = either usage (doProfileDistribution bam) $ argsToConf args
          handleOpt (_,    _,     []) = usage "Specify just one sorted, indexed BAM file"
          usage errs = do prog <- getProgName
                          let progline = prog ++ " [OPTIONS] <BAM>"
                          hPutStr stderr $ usageInfo progline optDescrs
                          hPutStrLn stderr errs

data TrxSeq = NoTrxSeq
            | TrxSeq !Int !BS.ByteString
              
trxSeq :: Maybe BS.ByteString -> Loc.ContigLoc -> TrxSeq
trxSeq Nothing _cdsloc = NoTrxSeq
trxSeq (Just sequ) cdsloc = case Loc.startPos cdsloc of
  (Pos.Pos off Plus) -> TrxSeq (fromIntegral off) sequ
  (Pos.Pos _off Minus) -> NoTrxSeq

trxCdsSeq :: TrxSeq -> (Int, Int) -> BS.ByteString
trxCdsSeq NoTrxSeq (cdsStart, cdsEnd) = BS.replicate (1 + cdsEnd - cdsStart) 'N'
trxCdsSeq (TrxSeq cdsOrigin sequ) (cdsStart, cdsEnd) = Loc.seqDataPad sequ trxLoc
  where trxLoc = Loc.fromStartEnd (fromIntegral $ cdsStart + cdsOrigin) (fromIntegral $ cdsEnd + cdsOrigin)

doProfileDistribution :: FilePath -> Conf -> IO ()
doProfileDistribution bamfile conf =
  withFile (concat [ confOut conf, "-pauses.txt" ]) WriteMode $ \hpause ->
  withFile (concat [ confOut conf, "-pausestat.txt" ]) WriteMode $ \hstat ->
  BamIndex.withIndex bamfile $ \bidx ->
  mapOverTranscripts ([ confBed conf ]) $ \trx -> 
  flip (maybe (return ())) (cds trx) $ \cdsloc -> do
    msequ <- maybe (return Nothing) (flip FaIdx.readLoc (location trx)) $! confFasta conf
    asites <- readASiteDelta $ confASite conf
    ntprof <- transcriptNtProfile asites bidx trx
    let cprof = profileFromStart Nothing cdsloc ntprof
        (statline, pauselines) = pauses conf trx (trxSeq msequ cdsloc) cprof
    hPutStrLn hstat statline
    mapM_ (hPutStrLn hpause) pauselines
    hFlush hpause

pauses :: Conf -> Transcript -> TrxSeq -> U.Vector Int -> (String, [String])
pauses conf trx tsequ p = ( intercalate "\t" $ kgid : fields, pauselines )
    where ( fields, pausefields ) = pauseFields conf tsequ p
          pauselines = map (intercalate "\t" . (kgid :)) pausefields
          kgid = BS.unpack . unSeqLabel . geneId $ trx

pauseFields :: Conf -> TrxSeq -> U.Vector Int -> ( [String], [[String]] )
pauseFields conf tsequ cdsvec | mostlyZero trimvec = ( [ "0" ], [] )
                              | median < minMedian = ( [ showFFloat (Just 1) median "" ], [] )
                              | iqr > maxIQR * median = ( [ showFFloat (Just 1) median ""
                                                          , showFFloat (Just 1) iqr ""
                                                          ]
                                                        , [] )
                              | otherwise = ( [ showFFloat (Just 1) median ""
                                              , showFFloat (Just 1) iqr ""
                                              , show . length $ pauselines
                                              ]
                                            , pauselines
                                            ) 
    where trimvec = U.slice 15 (U.length cdsvec - 20) cdsvec
          trimxvec = U.map fromIntegral trimvec
          median = continuousBy medianUnbiased 1 2 trimxvec
          iqr = midspread medianUnbiased 4 trimxvec
          minMedian = confMinMedian conf
          maxIQR = 3.0
          minPause = confMinPause conf
          stop = U.length cdsvec - 1
          pauselines = mapMaybe pauseline . U.toList . U.indexed $ cdsvec
          pauseline (idx, ct) 
              | and [ idx >= 5, idx < (stop - 1), fromIntegral ct >= minPause * median ]
                  = Just [ show idx, show (stop - idx), show ct
                         , showFFloat (Just 1) (fromIntegral ct / median) "" 
                         , pepsequ idx
                         , BS.unpack $! ntsequ idx
                         ]
              | otherwise = Nothing
          ntsequ aaidx = trxCdsSeq tsequ (aaidx * 3 - 30, aaidx * 3 + 17)
          pepsequ aaidx = let sequ = trxCdsSeq tsequ (aaidx * 3 - 30, aaidx * 3 + 17)
                          in [ oneLetterM $ ntToAaAt sequ (i * 3) | i <- [0..14] ]
          mostlyZero v = nZero * 2 > U.length v
            where nZero = U.foldl' countZero 0 v
                  countZero ni 0 = succ ni
                  countZero ni _ = ni

data Arg = ArgOutput { unArgOutput :: !String }
         | ArgBed { unArgBed :: !String }
         | ArgASite { unArgASite :: !String }
         | ArgFasta { unArgFasta :: !String }
         | ArgMinMedian { unArgMinMedian :: !String }
         | ArgMinPause {unArgMinPause :: !String }
         deriving (Show, Read, Eq, Ord)

argOutput :: Arg -> Maybe String
argOutput (ArgOutput del) = Just del
argOutput _ = Nothing

argBed :: Arg -> Maybe String
argBed (ArgBed bed) = Just bed
argBed _ = Nothing

argASite :: Arg -> Maybe String
argASite (ArgASite aSite) = Just aSite
argASite _ = Nothing

argFasta :: Arg -> Maybe String
argFasta (ArgFasta fasta) = Just fasta
argFasta _ = Nothing

argMinMedian :: Arg -> Maybe String
argMinMedian (ArgMinMedian minmedian) = Just minmedian
argMinMedian _ = Nothing

argMinPause :: Arg -> Maybe String
argMinPause (ArgMinPause minpause) = Just minpause
argMinPause _ = Nothing

optDescrs :: [OptDescr Arg]
optDescrs = [ Option ['o'] ["output"]     (ReqArg ArgOutput "OUTFILE")   "Output filename"
            , Option ['b'] ["bed"]        (ReqArg ArgBed "BED")          "Bed filename"
            , Option ['a'] ["asite"]      (ReqArg ArgASite "ASITEFILE")  "A site offsets filename"
            , Option ['f'] ["fasta"]      (ReqArg ArgFasta "FASTA")      "Indexed fasta sequence file"
            , Option []    ["min-median"] (ReqArg ArgMinMedian "MEDIAN") "Minimum per-codon median reads"
            , Option []    ["min-pause"]  (ReqArg ArgMinPause "PAUSE")   "Minimum reads at a pause, relative to median"
            ]

argsToConf :: [Arg] -> Either String Conf
argsToConf = runReaderT conf
    where conf = Conf <$> 
                 findBed <*>
                 findOutput <*>
                 findASite <*>
                 findFasta <*>
                 findMinMedian <*>
                 findMinPause
          findOutput = ReaderT $ maybe (Left "No out base") return  . listToMaybe . mapMaybe argOutput
          findBed = ReaderT $ maybe (Left "No bedfile") return . listToMaybe . mapMaybe argBed
          findASite = ReaderT $ maybe (Left "No A sites") return . listToMaybe . mapMaybe argASite
          findFasta = ReaderT $ return . listToMaybe . mapMaybe argFasta
          findMinMedian = ReaderT $ maybe (return defaultMinMedian) parseDouble . listToMaybe . mapMaybe argMinMedian
          findMinPause = ReaderT $ maybe (return defaultMinPause) parseDouble . listToMaybe . mapMaybe argMinPause
          parseDouble = AP.parseOnly AP.double . BS.pack