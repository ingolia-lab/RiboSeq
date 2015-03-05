{-# LANGUAGE OverloadedStrings, ScopedTypeVariables, FlexibleInstances #-}

module Main
       where 

import Control.Applicative
import Control.Exception
import Control.Monad.Reader
import Control.Monad.Trans.Resource
import qualified Data.ByteString.Char8 as BS
import Data.List (intercalate, maximumBy)
import Data.Maybe
import Data.Ord
import Foreign.Marshal.Utils
import Numeric
import System.Console.CmdTheLine
import System.IO

import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import qualified Data.Conduit as C
import qualified Data.Conduit.List as C

import qualified Bio.SamTools.Bam as Bam
import qualified Bio.SamTools.Conduit as Bam

import Bio.SeqLoc.Bed
import qualified Bio.SamTools.BamIndex as BamIndex
import qualified Bio.SeqLoc.LocMap as SLM
import Bio.SeqLoc.LocRepr
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.Position as Pos
import Bio.SeqLoc.Transcript

import Bio.RiboSeq.BamFile
import Bio.RiboSeq.Framing

doFpFraming :: Conf -> IO ()
doFpFraming conf = do
  trxmap <- readAnnotMap conf
  fsio <- fsioNew (confFlank conf) (confFlank conf) (confLengths conf) 
  Bam.withBamInFile (confBamInput conf) $ \hin ->
    withMaybeBamOutFile conf (Bam.inHeader hin) $ \mhout ->
    runResourceT $ C.runConduit $
    Bam.sourceHandle hin C.$$
    C.mapM_ (\bam -> let bamfr = bamFraming (confCdsBody conf) trxmap bam
                     in liftIO $ do
                       fsioIncr fsio bamfr (maybe (-1) fromIntegral $ Bam.queryLength bam)
                       case mhout of
                         Nothing -> return ()
                         Just hout -> do bannot <- Bam.addAuxZ bam tagFraming (framingAux bamfr)
                                         Bam.put1 hout bannot)
  fstats <- fsioFreeze fsio
  writeFile ((confOutput conf) ++ "_frame_length.txt") $ frameLenTable (fsBody fstats)
  writeFile ((confOutput conf) ++ "_around_start.txt") $ metagene2D (fsStart fstats)
  writeFile ((confOutput conf) ++ "_around_end.txt") $ metagene2D (fsEnd fstats)
  writeFile ((confOutput conf) ++ "_framing_stats.txt") $ framingStats fstats
  return ()
  
tagFraming :: String
tagFraming = "ZF"

framingAux (Left err) = show err
framingAux (Right (FpFraming start end frame gene))
  = intercalate "/" [ BS.unpack . unSeqLabel $ gene, showmz start, showmz end, showmz frame ]
  where showmz (Just z) = showSigned showInt 0 z ""
        showmz Nothing = "*"

readAnnotMap :: Conf -> IO (SLM.SeqLocMap Transcript)
readAnnotMap conf = do
  trxs <- concat <$> mapM readBedTranscripts (confBeds conf)
  let trxmap = SLM.locatableSeqLocMap defaultBinSize trxs
  hPutStrLn stderr $! "Read " ++ show (length trxs) ++ " annotations from " ++ show (confBeds conf)
  return $! trxmap
  where defaultBinSize = 100000

withMaybeBamOutFile :: Conf -> Bam.Header -> (Maybe Bam.OutHandle -> IO ()) -> IO ()
withMaybeBamOutFile conf inHeader f = case confAnnotate conf of
  Nothing -> f Nothing
  Just annotName -> Bam.withBamOutFile annotName inHeader $ \hout -> f (Just hout)

frameLenTable :: LengthFrame -> String
frameLenTable fr = unlines $ header : proflines
  where header = unfields [ "length", "fract", "N0", "N1", "N2", "p0", "p1", "p2", "info" ]
        proflines = map profline [0..(U.length (lfProfile fr V.! 0) - 1)]
          where total = V.sum . V.map U.sum . lfProfile $ fr
                profline l = let ln0 = (lfProfile fr V.! 0) U.! l
                                 ln1 = (lfProfile fr V.! 1) U.! l
                                 ln2 = (lfProfile fr V.! 2) U.! l
                                 lttl = ln0 + ln1 + ln2
                                 entropy = lengthFrameEntropy ( ln0, ln1, ln2 )
                                 info = logBase 2 3 - entropy
                             in unfields $ 
                                [ show $ l + lfMinLen fr
                                , showfract lttl total
                                , show ln0, show ln1, show ln2
                                , showfract ln0 lttl
                                , showfract ln1 lttl
                                , showfract ln2 lttl
                                , showFFloat (Just 2) info ""
                                ]
        showfract numer denom = showFFloat (Just 4) fract ""
          where fract :: Double
                fract = (fromIntegral numer) / (fromIntegral denom)

lengthFrameEntropy :: (Int, Int, Int) -> Double
lengthFrameEntropy (n0, n1, n2) = negate . sum $ map nEntropy [n0, n1, n2] 
  where nEntropy n = f * logBase 2 f
          where f = (fromIntegral n) / total
        total = fromIntegral $ n0 + n1 + n2
        
metagene2D :: LengthFrame -> String
metagene2D fr = unlines $ header : proflines
  where header = unfields $ [ "pos", "ttl" ] ++ [ show (i + lfMinLen fr) | i <- [0..(U.length (lfProfile fr V.! 0) - 1)] ]
        minpos = fromIntegral . lfMinPos $ fr
        proflines = map profline [0..((V.length . lfProfile $ fr) - 1)]
          where profline p = let vpos = lfProfile fr V.! p
                                 vttl = U.sum vpos
                             in unfields $ 
                                (show $ p + minpos) : 
                                map show ( vttl : U.toList vpos )

framingStats :: FramingStats -> String
framingStats fstats = unlines . map unfields $ stats
  where stats = [ [ "TOTAL", show . fsTotal $ fstats ]
                , [ "Start", show . V.sum . V.map U.sum . lfProfile . fsStart $ fstats ]
                , [ "Body",  show . V.sum . V.map U.sum . lfProfile . fsBody  $ fstats ]
                , [ "End",   show . V.sum . V.map U.sum . lfProfile . fsEnd   $ fstats ]
                ] ++
                [ [ show bf, show ( fsFailure fstats U.! (fromEnum bf) ) ] | (bf :: BamFailure) <- [minBound..maxBound] ]

unfields :: [String] -> String
unfields = intercalate "\t"

data Conf = Conf { confBamInput :: !FilePath
                 , confOutput :: !(FilePath) 
                 , confBeds :: ![FilePath]
                 , confFlank :: !(Pos.Offset, Pos.Offset)
                 , confCdsBody :: !(Pos.Offset, Pos.Offset)
                 , confLengths :: !(Int, Int)
                 , confAnnotate :: !(Maybe FilePath)
                 } deriving (Show)

argConf :: Term Conf
argConf = Conf <$>
          argBamInput <*>
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

argBamInput :: Term FilePath
argBamInput = required $ pos 0 Nothing $ posInfo
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
