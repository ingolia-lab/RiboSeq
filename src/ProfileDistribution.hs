{-# LANGUAGE TypeFamilies, FlexibleContexts, ScopedTypeVariables #-}

module Main where

import Control.Applicative
import Control.Monad
import Control.Monad.Error
import Control.Monad.Reader
import qualified Data.ByteString.Char8 as BS
import Data.List
import Data.Maybe
import Numeric
import System.Console.GetOpt
import System.Environment
import System.FilePath
import System.IO

import qualified Data.Attoparsec as AP
import qualified Data.Attoparsec.Char8 as AP
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U

import Statistics.Quantile
import Statistics.Sample

import qualified Bio.SamTools.BamIndex as BamIndex
import qualified Bio.SamTools.FaIdx as FaIdx
import Bio.SeqLoc.LocRepr
import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.Position as Pos
import Bio.SeqLoc.Strand
import Bio.SeqLoc.Transcript

import Bio.RiboSeq.BamFile
import Bio.RiboSeq.CodonAssignment
import Bio.RiboSeq.Counting

data Conf = Conf { confBed :: !FilePath
                 , confOut :: !FilePath
                 , confASite :: !FilePath
                 }
          deriving (Show)
        
main :: IO ()
main = getArgs >>= handleOpt . getOpt RequireOrder optDescrs
    where handleOpt (_,    _,         errs@(_:_)) = usage (unlines errs)
          handleOpt (args, [bam], []) = either usage (doProfileDistribution bam) $ argsToConf args
          handleOpt (_,    _,     []) = usage "Specify just one sorted, indexed BAM file"
          usage errs = do prog <- getProgName
                          let progline = prog ++ " [OPTIONS] <BAM>"
                          hPutStr stderr $ usageInfo progline optDescrs
                          hPutStrLn stderr errs

doProfileDistribution :: FilePath -> Conf -> IO ()
doProfileDistribution bamfile conf =
  withFile (concat [ confOut conf, "-pauses.txt" ]) WriteMode $ \hpause ->
  withFile (concat [ confOut conf, "-pausestat.txt" ]) WriteMode $ \hstat ->
  BamIndex.withIndex bamfile $ \bidx ->
  mapOverTranscripts ([ confBed conf ]) $ \trx -> 
  flip (maybe (return ())) (cds trx) $ \cdsloc -> do
    asites <- readASite $ confASite conf
    ntprof <- transcriptNtProfile (aSiteDelta asites) bidx trx
    let cprof = profileFromStart Nothing cdsloc ntprof
        (statline, pauselines) = pauses trx cprof
    hPutStrLn hstat statline
    mapM_ (hPutStrLn hpause) pauselines
    hFlush hpause

pauses :: Transcript -> U.Vector Int -> (String, [String])
pauses trx p = ( intercalate "\t" $ kgid : fields, pauselines )
    where ( fields, pausefields ) = pauseFields p
          pauselines = map (intercalate "\t" . (kgid :)) pausefields
          kgid = BS.unpack . unSeqLabel . geneId $ trx

pauseFields :: U.Vector Int -> ( [String], [[String]] )
pauseFields cdsvec | mostlyZero trimvec = ( [ "0" ], [] )
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
          minMedian = 3.0 - 1e6
          maxIQR = 3.0
          minPause = 25.0
          stop = U.length cdsvec - 1
          pauselines = mapMaybe pauseline . U.toList . U.indexed $ cdsvec
          pauseline (idx, ct) 
              | and [ idx >= 5, idx < (stop - 1), fromIntegral ct >= minPause * median ]
                  = Just [ show idx, show (stop - idx), show ct, showFFloat (Just 1) (fromIntegral ct / median) "" ]
              | otherwise = Nothing
          mostlyZero v = nZero * 2 > U.length v
            where nZero = U.foldl' countZero 0 v
                  countZero ni 0 = succ ni
                  countZero ni _ = ni

data Arg = ArgOutput { unArgOutput :: !String }
         | ArgBed { unArgBed :: !String }
         | ArgASite { unArgASite :: !String }
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

optDescrs :: [OptDescr Arg]
optDescrs = [ Option ['o'] ["output"]     (ReqArg ArgOutput "OUTFILE")   "Output filename"
            , Option ['b'] ["bed"]        (ReqArg ArgBed "BED")          "Bed filename"
            , Option ['a'] ["asite"]      (ReqArg ArgASite "ASITEFILE")  "A site offsets filename"
            ]

argsToConf :: [Arg] -> Either String Conf
argsToConf = runReaderT conf
    where conf = Conf <$> 
                 findBed <*>
                 findOutput <*>
                 findASite
          findOutput = ReaderT $ maybe (Left "No out base") return  . listToMaybe . mapMaybe argOutput
          findBed = ReaderT $ maybe (Left "No bedfile") return . listToMaybe . mapMaybe argBed
          findASite = ReaderT $ maybe (Left "No A sites") return . listToMaybe . mapMaybe argASite

