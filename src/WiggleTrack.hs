{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Main
       where 

import Control.Applicative
import Control.Exception
import Control.Monad.Reader
import qualified Data.ByteString.Char8 as BS
import Data.Int
import Data.List
import Data.Maybe
import Data.Ord
import Numeric
import System.Console.GetOpt
import System.Environment
import System.FilePath
import System.IO

import qualified Data.Attoparsec.Char8 as AP
import qualified Data.Conduit as C
import qualified Data.Conduit.List as C
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM

import qualified Bio.SamTools.Bam as Bam
import qualified Bio.SamTools.Conduit as Bam
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
    bracket (Bam.openBamInFile bam) Bam.closeInHandle $ \bin -> do
      aSiteDelta <- readASiteDelta $ confASite conf
      let asite = Bam.refSpLoc >=> aSitePos aSiteDelta
      countBam conf bin (hfwd, hrev) asite
                       
writeChrSizes :: [Bam.HeaderSeq] -> FilePath -> IO ()
writeChrSizes tseqs outname = withFile outname WriteMode $ \hout ->
  let writeChrSizeLine hseq = hPutStrLn hout $ concat [ BS.unpack . Bam.name $ hseq
                                                      , "\t"
                                                      , show . Bam.len $ hseq
                                                      ]
  in mapM_ writeChrSizeLine tseqs

data CountWindow = CountWindow { cwFwd :: !(UM.IOVector Int)
                               , cwRev :: !(UM.IOVector Int)
                               }
cwLength :: CountWindow -> Int
cwLength = UM.length . cwFwd

cwNew :: Int -> IO CountWindow
cwNew winlen = CountWindow <$> UM.replicate winlen 0 <*> UM.replicate winlen 0

cwCountOne :: CountWindow -> Int -> Strand -> IO ()
cwCountOne (CountWindow ctfwd ctrev) off strand = incr ctstrand off
  where ctstrand = case strand of 
          Plus -> ctfwd
          Minus -> ctrev
        incr v i = UM.read v i >>= UM.write v i . (succ $!)

data WindowCounts = WindowCounts { wcFwd, wcRev :: !(U.Vector Int)
                                 , wcOffset :: !Pos.Offset
                                 }

data CountWindowSet = CountWindowSet { cwsBefore, cwsCurr :: !CountWindow
                                     , cwsOffset :: !Pos.Offset
                                     }

cwsCountOne :: CountWindowSet -> Pos.Pos -> IO ()
cwsCountOne (CountWindowSet before curr winoff) (Pos.Pos off strand)
  = case fromIntegral $ off - winoff of
  x | x < negate (cwLength before) -> hPutStrLn stderr "cwsCountOne: skipping before"
    | x < 0 -> cwCountOne before (x + cwLength before) strand
    | x < cwLength curr -> cwCountOne curr x strand
    | otherwise -> hPutStrLn stderr "cwsCountOne: skipping after"

cwsNew :: Int -> Pos.Offset -> IO CountWindowSet
cwsNew winlen off0 = CountWindowSet <$>
                     cwNew winlen <*>
                     cwNew winlen <*>
                     pure off0

cwsBounds :: CountWindowSet -> (Pos.Offset, Pos.Offset, Pos.Offset)
cwsBounds (CountWindowSet before curr winoff)
  = let !start = (fromIntegral winoff) - (fromIntegral $ cwLength before)
        !end = start + fromIntegral (cwLength curr) - 1
        !next = start + (2 * fromIntegral (cwLength curr)) - 1
    in (start, end, next)

cwsAdvance :: CountWindowSet -> IO (CountWindowSet, WindowCounts)
cwsAdvance (CountWindowSet before@(CountWindow beforeFwd beforeRev) curr winoff) = do
  prevFwd <- U.freeze beforeFwd
  prevRev <- U.freeze beforeRev
  let !prev = WindowCounts prevFwd prevRev (winoff - fromIntegral (cwLength before))
  next <- cwNew (cwLength curr)
  let !advanced = CountWindowSet curr next (winoff + fromIntegral (cwLength curr))
  return (advanced, prev)

cwsFinish :: CountWindowSet -> IO WindowCounts
cwsFinish (CountWindowSet before@(CountWindow beforeFwd beforeRev) (CountWindow currFwd currRev) winoff) = do
  fstFwd <- U.freeze beforeFwd
  sndFwd <- U.freeze currFwd
  fstRev <- U.freeze beforeRev
  sndRev <- U.freeze currRev
  return $! WindowCounts (fstFwd U.++ sndFwd) (fstRev U.++ sndRev) (winoff - fromIntegral (cwLength before))

countBam :: Conf -> Bam.InHandle -> (Handle, Handle) -> (Bam.Bam1 -> Maybe Pos.Pos) -> IO ()
countBam conf bin hs asite = do cws0 <- cwsNew (confWindowSize conf) 0
                                (cws', _) <- Bam.sourceHandle bin C.$$ C.foldM countOne (cws0, -1)
                                cwsFinish cws' >>= wcPutWiggle hs conf
  where countOne (cws, tidx0) b = case Bam.targetID b of
          Nothing -> return (cws, tidx0)
          (Just tidx) | tidx < tidx0 -> error "countTarget: Decreasing target index, BAM unsorted?"
                      | tidx == tidx0 -> let (start, end, next) = cwsBounds cws
                                         in case asite b of
                                           Nothing -> return (cws, tidx0)
                                           (Just p) | Pos.offset p < start -> error "countTarget: Decreasing position, BAM unsorted?"
                                                    | Pos.offset p < end -> do cwsCountOne cws p
                                                                               return (cws, tidx)
                                                    | Pos.offset p < next -> do (cws', wc) <- cwsAdvance cws
                                                                                wcPutWiggle hs conf wc
                                                                                cwsCountOne cws' p
                                                                                return (cws', tidx)
                                                    | otherwise -> do cwsFinish cws >>= wcPutWiggle hs conf
                                                                      cws' <- cwsNew (confWindowSize conf) (Pos.offset p)
                                                                      cwsCountOne cws' p
                                                                      return (cws', tidx)
                      | otherwise -> do cwsFinish cws >>= wcPutWiggle hs conf
                                        let !newchr = Bam.targetSeqName (Bam.inHeader bin) tidx
                                        forM_ [fst hs, snd hs] $ \h -> hPutWiggleChrom h newchr
                                        cws' <- cwsNew (confWindowSize conf) (maybe 0 fromIntegral $ Bam.position b)
                                        countOne (cws', tidx) b -- N.B. re-enter with new target index

hPutWiggleChrom :: Handle -> BS.ByteString -> IO ()
hPutWiggleChrom h name = hPutStrLn h . unwords $ [ "variableStep", "chrom=" ++ BS.unpack name, "span=1" ]

wcPutWiggle :: (Handle, Handle) -> Conf -> WindowCounts -> IO ()
wcPutWiggle (hfwd, hrev) conf (WindowCounts fwd rev off) = do
  hPutWiggleWindow hfwd conf fwd off
  hPutWiggleWindow hrev conf rev off

hPutWiggleWindow :: Handle -> Conf -> U.Vector Int -> Pos.Offset -> IO ()
hPutWiggleWindow h conf ctvec off = U.imapM_ putDatum ctvec
  where putDatum i ct = case (fromIntegral off) + (fromIntegral i) of
          (x :: Int64) | x >= 0 && ct > 0 -> let qct = confQNorm conf * fromIntegral ct
                                             in hPutStrLn h $ shows (x + 1) . (' ' :) . showFFloat (Just 2) qct $ ""
                       | otherwise -> return ()

data Conf = Conf { confOutput :: !(FilePath) 
                 , confASite :: !FilePath
                 , confQNorm :: !Double
                 , confChrSizes :: !(Maybe FilePath)
                 , confWindowSize :: !Int
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
            , Option ['q'] ["qnorm"]      (ReqArg ArgQNorm "QNORM")       "Multiplicative scaling factor"
            , Option ['x'] ["chrsizes"]   (ReqArg ArgChrSizes "SIZEFILE") "Chromosome size output file"
            ]

argsToConf :: [Arg] -> Either String Conf
argsToConf = runReaderT conf
    where conf = Conf <$> 
                 findOutput <*>
                 findASite <*>
                 findQNorm <*>
                 findChrSizes <*>
                 pure defaultWindowSize
          findOutput = ReaderT $ maybe (Left "No out base") return  . listToMaybe . mapMaybe argOutput
          findASite = ReaderT $ maybe (Left "No A sites") return . listToMaybe . mapMaybe argASite
          findQNorm = ReaderT $ maybe (return 1.0) parseDouble . listToMaybe . mapMaybe argQNorm
            where parseDouble = AP.parseOnly AP.double . BS.pack
          findChrSizes = ReaderT $ return . listToMaybe . mapMaybe argChrSizes
          defaultWindowSize = 1024
