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
import qualified Bio.SeqLoc.Location as Loc
import qualified Bio.SeqLoc.LocMap as SLM
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.Position as Pos
import qualified Bio.SeqLoc.SpliceLocation as SpLoc
import Bio.SeqLoc.Strand
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
    inputSource conf C.$= C.mapMaybeM (annotate conf trxmap) C.$$ outputSink conf

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

annotate :: (MonadIO m) => Conf -> SLM.SeqLocMap Transcript -> Bam.Bam1 -> m (Maybe Bam.Bam1)
annotate conf trxmap b = liftIO $ do
  b1 <- if null hits
        then return b
        else Bam.addAuxZ b trxNameTag $ annotNames conf hits
  b2 <- if null hits
        then return b1
        else Bam.addAuxZ b1 trxPosTag $ annotPoses conf hits
  b3 <- if null hits || not (cAnnotStrand conf)
        then return b2
        else Bam.addAuxZ b2 trxStrandTag $ annotStrands conf hits
  return $ Just b3
  where trxNameTag = "ZT"
        trxPosTag = "ZP"
        trxStrandTag = "ZS"
        hits = trxHits conf trxmap b

annotNames :: Conf -> [(Transcript, Loc.ContigLoc)] -> String
annotNames _conf hits = intercalate "," . map (BS.unpack . unSeqLabel . trxId . fst) $ hits

annotPoses :: Conf -> [(Transcript, Loc.ContigLoc)] -> String
annotPoses _conf hits = intercalate "," . map (show . Pos.unOff . Loc.offset5 . snd) $ hits

annotStrands :: Conf -> [(Transcript, Loc.ContigLoc)] -> String
annotStrands _conf hits = intercalate "," . map (strandchr . Loc.strand . snd) $ hits
  where strandchr str = case str of Plus -> "+"; Minus -> "-"

trxHits :: Conf -> SLM.SeqLocMap Transcript -> Bam.Bam1 -> [(Transcript, Loc.ContigLoc)]
trxHits conf trxmap b = case Bam.refSeqLoc b of
  Just bamloc -> mapMaybe toContig $ SLM.queryLocInto (cStrandSpecific conf) bamloc trxmap
  Nothing -> []
  where toContig (trx, sploc) = case Loc.toContigs sploc of
          [c1] -> Just (trx, c1)
          _ -> Nothing

data Conf = Conf { cBamInput :: !FilePath
                 , cBedInput :: !FilePath
                 , cBamOutput :: !FilePath
                 , cBinSize :: !Int
                 , cTamInput :: !Bool
                 , cTamOutput :: !Bool
                 , cAnnotStrand :: !Bool
                 , cStrandSpecific :: !(Maybe Strand)
                 } deriving (Show)

argConf :: Term Conf
argConf = Conf <$>
          argBamInput <*>
          argBedInput <*>
          argBamOutput <*>
          argBinSize <*>
          argTamInput <*>
          argTamOutput <*>
          argAnnotStrand <*> 
          argStrandSpecific

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

argAnnotStrand :: Term Bool
argAnnotStrand = value $ flag $ (optInfo ["a", "annotate-strand"]) { optDoc = "Annotate relative strand of feature overlap" }

argStrandSpecific :: Term (Maybe Strand)
argStrandSpecific = value $ vFlag Nothing [ ( Just Plus,  (optInfo ["f", "forward"]) { optDoc = "Annotate only \"forward\" (i.e., same strand) overlaps" })
                                          , ( Just Minus, (optInfo ["r", "reverse"]) { optDoc = "Annotate only \"reverse\" (i.e., opposite strand) overlaps" })
                                          ]

                    
