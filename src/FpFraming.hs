{-# LANGUAGE OverloadedStrings, ScopedTypeVariables, FlexibleInstances #-}

module Main
       where 

import Control.Applicative
import Control.Exception
import Control.Monad.Reader
import qualified Data.ByteString.Char8 as BS
import Data.List (intercalate, maximumBy)
import Data.Maybe
import Data.Ord
import Foreign.Marshal.Utils
import Numeric
import System.Console.CmdTheLine
import System.IO

import qualified Data.Vector as V

import qualified Bio.SamTools.BamIndex as BamIndex
import Bio.SeqLoc.LocRepr
import qualified Bio.SeqLoc.Position as Pos

import Bio.RiboSeq.BamFile
import Bio.RiboSeq.Framing

doFpFraming :: Conf -> IO ()
doFpFraming conf = do trstartio <- newTerminus (confFlank conf) (confLengths conf)
                      trendio <- newTerminus (confFlank conf) (confLengths conf)
                      frameio <- newFrame (confLengths conf)
                      let count trx iloc = do countAtStart trstartio trx iloc
                                              countAtEnd trendio trx iloc
                                              countInCds (confCdsBody conf) frameio trx iloc
                      withMany (\bam -> bracket (BamIndex.open bam) (BamIndex.close)) (confBamInputs conf) $ \bidxs ->
                        mapOverTranscripts (confBeds conf) $ \trx ->
                        forM_ bidxs $ \bidx -> mapOverBams bidx (count trx) trx
                      trstart <- freezeTerminus trstartio
                      trend <- freezeTerminus trendio
                      frame <- freezeFrame frameio
                      writeFile (confOutput conf ++ "_start_pos_len.txt") . posLenTable $ trstart
                      writeFile (confOutput conf ++ "_end_pos_len.txt") . posLenTable $ trend
                      writeFile (confOutput conf ++ "_frame_len.txt") . frameLenTable $ frame
                      writeFile (confOutput conf ++ "_asite_report.txt") $ framingTable frame trstart trend

frameLenTable :: Framing -> String
frameLenTable fr = unlines $ header ++ proflines (frprofile fr)
  where header = [ unwords [ "# lengths ", show . frminlen $ fr, show . frmaxlen $ fr ]
                 ]
        proflines prof = V.toList . V.imap profline $ prof
          where total = V.sum . V.map V.sum $ prof
                profline i1 v = let ltotal = V.sum v
                                    framect f = show $ v V.! f
                                    framefract f = showfract (v V.! f) ltotal
                                in unfields $ 
                                   [ show $ i1 + frminlen fr
                                   , showfract ltotal total
                                   ] 
                                   ++ map framect [0..2]
                                   ++ map framefract [0..2]
        showfract numer denom = showFFloat (Just 4) fract ""
          where fract :: Double
                fract = (fromIntegral numer) / (fromIntegral denom)

minLenFract :: Double
minLenFract = 0.05

startRange :: (Pos.Offset, Pos.Offset)
startRange = (-17, -8)

startShift :: Pos.Offset
startShift = 3

endRange :: (Pos.Offset, Pos.Offset)
endRange = (-22, -13)

endShift :: Pos.Offset
endShift = -2

framingTable :: Framing -> Terminus -> Terminus -> String
framingTable fr trstart trend = unlines . (header : ) . map framingLine $ wantedLength
  where wantedLength = filter ((> minLenFract) . lenFract fr) [(frminlen fr)..(frmaxlen fr)]
        header = unfields [ "# Len", "Fract", "AtStart", "AtEnd", "Info", "Frame0", "Frame1", "Frame2", "ASite" ]
        framingLine len = unfields $ [ show len
                                     , showFFloat (Just 2) (lenFract fr len) ""
                                     , BS.unpack . repr $ startPeak
                                     , BS.unpack . repr $ endPeak
                                     , showFFloat (Just 2) frameInfo ""
                                     ] ++ map showfr [0..2] ++ [ bestASite ]
          where startPeak = startShift + (negate $ terminusPeak startRange trstart len)
                endPeak = endShift + (negate $ terminusPeak endRange trend len)
                nfr f = (frprofile fr V.! (len - frminlen fr)) V.! f
                frameInfo = logBase 2 3 - lengthFrameEntropy (nfr 0, nfr 1, nfr 2)
                showfr f = showFFloat (Just 2) (fromIntegral (nfr f) / lenttl) ""
                lenttl = (fromIntegral $ V.sum $ frprofile fr V.! (len - frminlen fr)) :: Double
                startfr = fromIntegral $ (negate startPeak) `mod` 3
                endfr = fromIntegral $ (negate endPeak) `mod` 3
                bestfr = maximumBy (comparing nfr) [0..2]
                bestASite | startPeak == endPeak && bestfr == startfr = BS.unpack . repr $ startPeak
                          | endfr == bestfr = BS.unpack . repr $  endPeak
                          | startfr == bestfr = BS.unpack . repr $ startPeak
                          | startPeak == endPeak = (BS.unpack . repr $ startPeak) ++ "???"
                          | otherwise = "***"

lengthFrameEntropy :: (Int, Int, Int) -> Double
lengthFrameEntropy (n0, n1, n2) = negate . sum $ map nEntropy [n0, n1, n2] 
  where nEntropy n = f * logBase 2 f
          where f = (fromIntegral n) / total
        total = fromIntegral $ n0 + n1 + n2
        
posLenTable :: Terminus -> String
posLenTable tr = unlines $ header ++ proflines (profile tr)
  where header = [ unwords [ "# positions ", show . Pos.unOff .  before $ tr, show . Pos.unOff . after $ tr]
                 , unwords [ "# lengths ", show . minlen $ tr, show . maxlen $ tr ]
                 ]        
        proflines = V.toList . V.imap profline
        b = fromIntegral . before $ tr
        profline i1 v = unfields $ (show $ i1 + b) : show ttl : (V.toList . V.map show $ v)
          where ttl = V.sum v


unfields :: [String] -> String
unfields = intercalate "\t"

data Conf = Conf { confBamInputs :: [FilePath]
                 , confOutput :: !(FilePath) 
                 , confBeds :: [FilePath]
                 , confFlank :: !(Pos.Offset, Pos.Offset)
                 , confCdsBody :: !(Pos.Offset, Pos.Offset)
                 , confLengths :: !(Int, Int)
                 , confAnnotate :: !(Maybe FilePath)
                 } deriving (Show)

argConf :: Term Conf
argConf = Conf <$>
          argBamInputs <*>
          argOutput <*>
          argBedFiles <*>
          argFlank <*>
          argBody <*>
          argLengths <*>
          argAnnotate

instance ArgVal Pos.Offset where
  converter = let (intParser :: ArgParser Int, intPrinter :: ArgPrinter Int) = converter
              in ( either Left (Right . fromIntegral) . intParser
                 , intPrinter . fromIntegral
                 )

instance ArgVal (Pos.Offset, Pos.Offset) where
  converter = pair ','

instance ArgVal (Int, Int) where
  converter = pair ','

argBamInputs :: Term [FilePath]
argBamInputs = nonEmpty $ posAny [] $ posInfo
  { posName = "BAM", posDoc = "BAM format alignment file" }

argOutput :: Term FilePath
argOutput = required $ opt Nothing $ (optInfo ["o", "output"])
  { optName = "OUTBASE", optDoc = "Base filename for output files" }

argBedFiles :: Term [FilePath]
argBedFiles = nonEmpty $ optAll [] $ (optInfo ["b", "bed"])
  { optName = "BED", optDoc = "BED-format annotation filename" }

argFlank :: Term (Pos.Offset, Pos.Offset)
argFlank = value $ opt (-100, 100) $ (optInfo ["f", "flanking"])
  { optName = "START,END", optDoc = "Range of profiles surrounding the start and end codons" }

argBody :: Term (Pos.Offset, Pos.Offset)
argBody = value $ opt (34, 31) $ (optInfo ["c", "cdsbody"])
  { optName = "AFTERSTART,BEFOREEND", optDoc = "Offsets from the start and end of the gene for framing analysis" }

argLengths :: Term (Int, Int)
argLengths = value $ opt (26, 34) $ (optInfo ["l", "lengths"])
  { optName = "MINLEN,MAXLEN", optDoc = "Length frange for framing analysis" }

argAnnotate :: Term (Maybe FilePath)
argAnnotate = value $ opt Nothing $ (optInfo ["a", "annotate"])
  { optName = "ANNOTATED.BAM", optDoc = "Write output BAM file annotated wiht framing information" }

main :: IO ()
main = run ( fpframe, info )
  where fpframe = doFpFraming <$> argConf
        info = defTI { termName = "fp-framing"
                     , version = "150304"
                     , termDoc = "Calculates ribosome profiling QC information -- reading frame bias and start and stop codon meta-genes"
                     }
