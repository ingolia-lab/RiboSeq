{-# LANGUAGE OverloadedStrings #-}

module Main
       where 

import Control.Applicative
import Control.Exception
import Control.Monad.Reader
import qualified Data.ByteString.Char8 as BS
import Data.List
import Data.Maybe
import Data.Ord
import Numeric
import System.Console.GetOpt
import System.Environment
import System.FilePath
import System.IO

import qualified Data.Attoparsec.Char8 as AP
import qualified Data.Iteratee as Iter
import qualified Data.Vector.Unboxed.Mutable as UM

import qualified Bio.SamTools.Bam as Bam
import qualified Bio.SamTools.BamIndex as BamIndex
import qualified Bio.SamTools.Iteratee as BamIter
import qualified Bio.SeqLoc.Location as Loc
import qualified Bio.SeqLoc.Position as Pos
import Bio.SeqLoc.Strand

import Bio.RiboSeq.CodonAssignment

main :: IO ()
main = getArgs >>= handleOpt . getOpt RequireOrder optDescrs
    where handleOpt (_,    _,       errs@(_:_)) = usage (unlines errs)
          handleOpt (_,    [],      [])         = usage "Specify a BAM file"
          handleOpt (args, [bam],   [])         = either usage (doWiggleTrack bam) $ argsToConf args
          handleOpt (_,    (_:_:_), [])         = usage "Specify just one BAM file"
          usage errs = do prog <- getProgName
                          let progline = prog ++ " [OPTIONS] <BAM>"
                          hPutStr stderr $ usageInfo progline optDescrs
                          hPutStrLn stderr errs

doWiggleTrack :: FilePath -> Conf -> IO ()
doWiggleTrack bam conf = do 
  tseqs <- Bam.withBamInFile bam $ return . Bam.targetSeqList . Bam.inHeader
  maybe (return ()) (writeChrSizes tseqs) $ confChrSizes conf
  withFile (confOutputPlus conf) WriteMode $ \hfwd ->
    withFile (confOutputRev conf) WriteMode $ \hrev ->
    forM_ (sortBy (comparing Bam.len) tseqs) $ \tseq -> do
      ct@(Count name ctfwd ctrev) <- targetSeqCounts conf tseq
      bracket (BamIndex.open bam) BamIndex.close $ \bidx -> do
        verbose conf $ unwords [ "    counting ", show bam ]
        countBam conf bidx ct
        verbose conf $ unwords [ "    done counting ", show bam ]
      verbose conf $ unwords ["  writing ", show $ confOutputPlus conf ]
      hPutWiggleChr hfwd conf name ctfwd
      verbose conf $ unwords ["  writing ", show $ confOutputRev conf ]      
      hPutWiggleChr hrev conf name ctrev
      verbose conf $ unwords ["  done writing ", show $ Bam.name tseq ]
                       
writeChrSizes :: [Bam.HeaderSeq] -> FilePath -> IO ()
writeChrSizes tseqs outname = withFile outname WriteMode $ \hout ->
  let writeChrSizeLine hseq = hPutStrLn hout $ concat [ BS.unpack . Bam.name $ hseq
                                                      , "\t"
                                                      , show . Bam.len $ hseq
                                                      ]
  in mapM_ writeChrSizeLine tseqs

data Count = Count { ctName :: !BS.ByteString, ctFwd :: !(UM.IOVector Int), ctRev:: !(UM.IOVector Int) }
data CountWindow = CountWindow { cwFwd :: !(UM.IOVector Int)
                               , cwRev :: !(UM.IOVector Int)
                               }
cwLength :: CountWindow -> Int
cwLength = UM.length . cwFwd

cwCountOne :: CountWindow -> Int -> Strand -> IO ()
cwCountOne (CountWindow ctfwd ctrev) off strand = incr ctstrand off
  where ctstrand = case strand of 
          Plus -> ctfwd
          Minus -> ctrev
        incr v i = UM.read v i >>= UM.write v i . (succ $!)

data CountWindowSet = CountWindowSet { cwsBefore, cwsWindow, cwsAfter :: !CountWindow
                                     , cwsOffset :: !Pos.Offset
                                     }

cwsCountOne :: CountWindowSet -> Pos.Pos -> IO ()
cwsCountOne (CountWindowSet before curr after winoff) (Pos.Pos off strand)
  = case fromIntegral $ off - winoff of
  x | x < negate (cwLength before) -> hPutStrLn stderr "cwsCountOne: skipping before"
    | x < 0 -> cwCountOne before (x + cwLength before) strand
    | x < cwLength curr -> cwCountOne curr x strand
    | x < (cwLength curr + cwLength after) -> cwCountOne after (x - cwLength curr) strand
    | otherwise -> hPutStrLn stderr "cwsCountOne: skipping after"

ctLength :: Count -> Int
ctLength = UM.length . ctFwd             

countPos :: Count -> Pos.Pos -> IO ()
countPos (Count _name ctfwd ctrev) (Pos.Pos off strand) = incr ctstrand $ fromIntegral off
  where ctstrand = case strand of 
          Plus -> ctfwd
          Minus -> ctrev
        incr v i = UM.read v i >>= UM.write v i . (succ $!)

targetSeqCounts :: Conf -> Bam.HeaderSeq -> IO Count
targetSeqCounts conf hseq = do verbose conf $ "Preparing to initialize count vectors"
                               verbose conf . unwords $ [ "  forward", BS.unpack $ Bam.name hseq
                                                        , show (Bam.len hseq), "positions"
                                                        ]
                               ctfwd <- newCountVector (Bam.len hseq)
                               verbose conf . unwords $ [ "  reverse", BS.unpack $ Bam.name hseq
                                                        , show (Bam.len hseq), "positions"
                                                        ]
                               ctrev <- newCountVector (Bam.len hseq)
                               verbose conf . unwords $ [ "Initialized", BS.unpack $ Bam.name hseq
                                                        , show (Bam.len hseq), "positions"
                                                        ]
                               return $! Count (Bam.name hseq) ctfwd ctrev
  where newCountVector len = UM.replicate (fromIntegral len) 0

countBam :: Conf -> BamIndex.IdxHandle -> Count -> IO ()
countBam conf bidx ct = do
  let noTarget = "Could not find target " ++ show (ctName ct)
  tidx <- maybe (error noTarget) return $! Bam.lookupTarget (BamIndex.idxHeader bidx) (ctName ct) 
  docount <- maybe (return countCoverage) countASite $ confASite conf
  BamIter.enumIndexRegion bidx tidx (0, fromIntegral $ ctLength ct - 1) (Iter.mapM_ docount) >>= Iter.run
  where countASite asitefile = readASiteDelta asitefile >>= \asites ->
          return $ maybe (return ()) (countPos ct) . (Bam.refSpLoc >=> aSitePos asites)
        countCoverage = maybe (return ()) countAll . Bam.refSpLoc
        countAll = mapM_ (countPos ct) . Loc.allPos

hPutWiggleChr :: Handle -> Conf -> BS.ByteString -> UM.IOVector Int -> IO ()
hPutWiggleChr h conf name ctvec = putHeader >> putData
  where putHeader = hPutStrLn h . unwords $ [ "variableStep", "chrom=" ++ BS.unpack name, "span=1" ]
        putData = forM_ [0..(UM.length ctvec - 1)] $ \idx -> 
          do ct <- UM.read ctvec idx
             let qct = confQNorm conf * fromIntegral ct
             when (ct > 0) $ hPutStrLn h $ shows (idx + 1) . (' ' :) . showFFloat (Just 2) qct $ ""

data Conf = Conf { confOutput :: !(FilePath) 
                 , confASite :: !(Maybe FilePath)
                 , confQNorm :: !Double
                 , confChrSizes :: !(Maybe FilePath)
                 } deriving (Show)

confOutputPlus :: Conf -> FilePath
confOutputPlus conf = case splitExtension $ confOutput conf of
  ( base, ext ) -> (base ++ "_fwd") `addExtension` ext
  
confOutputRev :: Conf -> FilePath
confOutputRev conf = case splitExtension $ confOutput conf of
  ( base, ext ) -> (base ++ "_rev") `addExtension` ext

verbose :: Conf -> String -> IO ()
verbose _conf message = hPutStrLn stderr message

data Arg = ArgOutput { unArgOutput :: !String }
         | ArgASite { unArgASite :: !String }
         | ArgQNorm { unArgQNorm :: !String }
         | ArgCoverage
         | ArgChrSizes { unArgChrSizes :: !String }
         deriving (Show, Read, Eq, Ord)

argOutput :: Arg -> Maybe String
argOutput (ArgOutput out) = Just out
argOutput _ = Nothing

argASite :: Arg -> Maybe String
argASite (ArgASite aSite) = Just aSite
argASite _ = Nothing

argQNorm :: Arg -> Maybe String
argQNorm (ArgQNorm qNorm) = Just qNorm
argQNorm _ = Nothing

argChrSizes :: Arg -> Maybe String
argChrSizes (ArgChrSizes chrSize) = Just chrSize
argChrSizes _ = Nothing

optDescrs :: [OptDescr Arg]
optDescrs = [ Option ['o'] ["output"]     (ReqArg ArgOutput "OUTFILE")    "Output filename"
            , Option ['a'] ["asite"]      (ReqArg ArgASite "ASITEFILE")   "A site offsets filename"
            , Option ['c'] ["coverage"]   (NoArg ArgCoverage)             "Total read coverage"
            , Option ['q'] ["qnorm"]      (ReqArg ArgQNorm "QNORM")       "Multiplicative scaling factor"
            , Option ['x'] ["chrsizes"]   (ReqArg ArgChrSizes "SIZEFILE") "Chromosome size output file"
            ]

argsToConf :: [Arg] -> Either String Conf
argsToConf = runReaderT conf
    where conf = Conf <$> 
                 findOutput <*>
                 findASite <*>
                 findQNorm <*>
                 findChrSizes
          findOutput = ReaderT $ maybe (Left "No out base") return  . listToMaybe . mapMaybe argOutput
          findASite = ReaderT $ \args ->
            let masites = listToMaybe . mapMaybe argASite $ args
            in if ArgCoverage `elem` args
               then maybe (return Nothing) (const $ Left "Cannot combine A sites and coverage") masites
               else maybe (Left "No A sites specified (and not coverage mode)") (return . Just) masites
          findQNorm = ReaderT $ maybe (return 1.0) parseDouble . listToMaybe . mapMaybe argQNorm
            where parseDouble = AP.parseOnly AP.double . BS.pack
          findChrSizes = ReaderT $ return . listToMaybe . mapMaybe argChrSizes
