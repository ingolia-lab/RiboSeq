{-# LANGUAGE BangPatterns, OverloadedStrings, RankNTypes #-}
module Main
       where

import Control.Applicative
import Control.Monad
import Control.Monad.Trans.Resource
import Control.Monad.Reader
import Control.Monad.Error
import qualified Data.ByteString.Char8 as BS
import qualified Data.Conduit as C
import qualified Data.Conduit.List as C
import qualified Data.HashSet as HS
import Data.List
import Data.Maybe
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import Numeric
import System.Console.CmdTheLine
import System.Environment
import System.IO

import Bio.SamTools.Bam as Bam
import Bio.SamTools.Conduit as Bam
import Bio.SeqLoc.Bed
import qualified Bio.SeqLoc.LocMap as SLM
import Bio.SeqLoc.Transcript

main :: IO ()
main = run ( bamannot, info )
  where info = defTI { termName = "bam-annotate"
                     , version = "0.0"
                     , termDoc = "Annotate BAM file with overlap of BED features"
                     }
        bamannot = bamAnnotate <$> argConf

bamAnnotate :: Conf -> IO ()
bamAnnotate conf = do
  trxmap <- readAnnotMap conf
  runResourceT $ C.runConduit $
    inputSource conf C.$= C.mapMaybe (annotate conf trxmap) C.$$ outputSink conf

readAnnotMap :: Conf -> IO (SLM.SeqLocMap Transcript)
readAnnotMap conf = do
  trxs <- readBedTranscripts (cBedInput conf)
  let trxmap = SLM.locatableSeqLocMap (fromIntegral . cBinSize $ conf) trxs
  hPutStrLn stderr $! "Read " ++ show (length trxs) ++ " annotations from " ++ show (cBedInput conf)
  return $! trxmap

inputSource :: (MonadResource m) => Conf -> C.Producer m Bam.Bam1
inputSource conf | cTamInput conf = Bam.sourceTamInFile (cBamInput conf)
                 | otherwise      = Bam.sourceBamInFile (cBamInput conf)

outputSink :: (MonadResource m) => Conf -> C.Consumer Bam.Bam1 m ()
outputSink conf | cTamOutput conf = Bam.sinkTamOutFile (cBamOutput conf)
                | otherwise       = Bam.sinkBamOutFile (cBamOutput conf)

annotate :: Conf -> SLM.SeqLocMap Transcript -> Bam.Bam1 -> Maybe Bam.Bam1
annotate _conf _trxmap b = Just b

data Conf = Conf { cBamInput :: !FilePath
                 , cBedInput :: !FilePath
                 , cBamOutput :: !FilePath
                 , cBinSize :: !Int
                 , cTamInput :: !Bool
                 , cTamOutput :: !Bool
                 } deriving (Show)

argConf :: Term Conf
argConf = Conf <$>
          argBamInput <*>
          argBedInput <*>
          argBamOutput <*>
          argBinSize <*>
          argTamInput <*>
          argTamOutput

argBamInput :: Term FilePath
argBamInput = required $ opt Nothing $ (optInfo ["i", "input"])
  { optName = "INPUT.BAM", optDoc = "BAM format input file" }

argBedInput :: Term FilePath
argBedInput = required $ opt Nothing $ (optInfo ["b", "bed"])
  { optName = "ANNOTATE.BED", optDoc = "BED format annotation file" }

argBamOutput :: Term FilePath
argBamOutput = required $ opt Nothing $ (optInfo ["o", "output"])
  { optName = "OUTPUT.BAM", optDoc = "BAM format output file" }

argBinSize :: Term Int
argBinSize = required $ opt defaultBinSize $ (optInfo ["z", "binsize"])
  { optName = "BIN-SIZE", optDoc = "Bin size for annotation map" }
  where defaultBinSize = Just 100000

argTamInput :: Term Bool
argTamInput = value $ flag $ (optInfo ["u", "text-input"]) { optDoc = "Text format input" }

argTamOutput :: Term Bool
argTamOutput = value $ flag $ (optInfo ["t", "text-output"]) { optDoc = "Text format output" }
