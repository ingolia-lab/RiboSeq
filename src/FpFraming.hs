{-# LANGUAGE OverloadedStrings, ScopedTypeVariables, FlexibleInstances #-}

module Main
       where 

import Control.Applicative
import Control.Monad.Reader
import Control.Monad.Trans.Resource
import qualified Data.ByteString.Char8 as BS
import Data.List (intercalate)
import Data.Maybe
import Numeric
import System.Console.CmdTheLine
import System.IO

import qualified Data.HashMap.Strict as HM
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import qualified Data.Conduit as C
import qualified Data.Conduit.List as C

import qualified Bio.SamTools.Bam as Bam
import qualified Bio.SamTools.Conduit as Bam

import Bio.SeqLoc.Bed
import qualified Bio.SeqLoc.LocMap as SLM
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.Position as Pos
import Bio.SeqLoc.Transcript

import Bio.RiboSeq.Framing

doFpFraming :: Conf -> IO ()
doFpFraming conf = do
  trxmap <- readAnnotMap conf
  fsio <- fsioNew (confFlank conf) (confFlank conf) (confLengths conf) 
  Bam.withBamInFile (confBamInput conf) $ \hin ->
    withMaybeBamOutFile conf (Bam.inHeader hin) $ \mhout ->
    runResourceT $ C.runConduit $
    Bam.sourceHandle hin C.$$
    C.mapM_ (\bam -> let bamfr = case liftM fromIntegral $ Bam.queryLength bam of
                           Nothing -> Left $ BamNoHit
                           Just len | len < confMinLength conf -> Left $ BamTooShort
                                    | len > confMaxLength conf -> Left $ BamTooLong
                                    | otherwise -> bamFraming (confCdsBody conf) trxmap bam
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

framingAux :: BamFramingResult -> String
framingAux (Left (BamFpFailure bf)) = show bf
framingAux (Left bf) = show bf
framingAux (Right (FpFraming start end frame gene))
  = intercalate "/" [ BS.unpack . unSeqLabel $ gene, showmz start, showmz end, showmz frame ]
  where showmz (Just z) = showSigned showInt 0 z ""
        showmz Nothing = "*"

readAnnotMap :: Conf -> IO (SLM.SeqLocMap Transcript)
readAnnotMap conf = do
  trxsRaw <- concat <$> mapM readBedTranscripts (confBedFiles conf)
  hPutStrLn stderr $! "Read " ++ show (length trxsRaw) ++ " annotations from " ++ show (confBedFiles conf)
  trxToGene <- mapM BS.readFile (confGeneFiles conf) >>=
               parseGeneMap . concat . map BS.lines
  hPutStrLn stderr $! "Read " ++ show (HM.size trxToGene) ++ " transcript-to-gene mappings from " ++ show (confGeneFiles conf)
  let trxs = mapMaybe (remapGenes trxToGene) trxsRaw
      ngene = HM.size . HM.fromList . map (\t -> (geneId t, ())) $ trxs
  hPutStrLn stderr $! "After remapping " ++ (show ngene) ++ " distinct genes in " ++ (show $ length trxs) ++ " transcripts"
  let trxmap = SLM.locatableSeqLocMap (confBinSize conf) trxs
  return $! trxmap
  where parseGeneMap ls = (HM.fromList . catMaybes) <$> mapM parseGeneMapping ls
        parseGeneMapping l
          | BS.null l = return Nothing
          | BS.isPrefixOf "#" l = return Nothing
          | otherwise = case BS.split '\t' l of
            [trxid,geneid] -> return . Just $! (BS.copy trxid, Just $ BS.copy geneid)
            [trxid] -> return . Just $! (BS.copy trxid, Nothing)
            _ -> do hPutStrLn stderr $ "Malformed transcript-to-gene mapping " ++ show l
                    return Nothing
        remapGenes trxToGene t = case HM.lookup (unSeqLabel . trxId $ t) trxToGene of
          Nothing -> Just t
          Just Nothing -> Nothing
          Just (Just g) -> Just $! t { geneId = toSeqLabel g }
        
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
  where stats = concat
                [ [ [ "TOTAL", "", show . fsTotal $ fstats ] ]
                , [ [ "", show bf,
                      show (fsFailure fstats U.! (fromEnum bf)),
                      showfract  (fsFailure fstats U.! (fromEnum bf)) (fsTotal fstats) ]
                  | bf <- badAlignment ]
                , [ [ "BadAlignment",  "", show badAlignCount, showfract badAlignCount (fsTotal fstats) ]
                  , [ "GoodAlignment", "", show goodAlignCount, showfract goodAlignCount (fsTotal fstats) ] ]
                , [ [ "", drop 2 $ show ff,
                      show (fsFailure fstats U.! (fromEnum $ BamFpFailure ff)),
                      showfract (fsFailure fstats U.! (fromEnum $ BamFpFailure ff)) (fsTotal fstats),
                      showfract (fsFailure fstats U.! (fromEnum $ BamFpFailure ff)) goodAlignCount ]
                  | ff <- [minBound..maxBound] ]
                , [ [ "BadAnnotation",  "", show badAnnotCount, showfract badAnnotCount (fsTotal fstats), showfract badAnnotCount goodAlignCount ]
                  , [ "GoodAnnotation", "", show goodAnnotCount, showfract goodAnnotCount (fsTotal fstats), showfract goodAnnotCount goodAlignCount ] ]                  
                , [ [],
                    [ "Start", show . V.sum . V.map U.sum . lfProfile . fsStart $ fstats ]
                  , [ "Body",  show . V.sum . V.map U.sum . lfProfile . fsBody  $ fstats ]
                  , [ "End",   show . V.sum . V.map U.sum . lfProfile . fsEnd   $ fstats ] ]
                ]
        badAlignCount = sum [ (fsFailure fstats) U.! (fromEnum bf) | bf <- badAlignment ]
        goodAlignCount = (fsTotal fstats) - badAlignCount
        badAnnotCount = sum [ (fsFailure fstats) U.! (fromEnum bf) | bf <- badAnnotation ]
        goodAnnotCount = goodAlignCount - badAnnotCount

showfract :: Int -> Int -> String
showfract numer denom = showFFloat (Just 4) fract ""
  where fract :: Double
        fract = (fromIntegral numer) / (fromIntegral denom)

unfields :: [String] -> String
unfields = intercalate "\t"

data Conf = Conf { confBamInput :: !FilePath
                 , confOutput :: !(FilePath) 
                 , confBedFiles :: ![FilePath]
                 , confGeneFiles :: ![FilePath]
                 , confFlank :: !(Pos.Offset, Pos.Offset)
                 , confCdsBody :: !(Pos.Offset, Pos.Offset)
                 , confLengths :: !(Int, Int)
                 , confAnnotate :: !(Maybe FilePath)
                 , confBinSize :: !Pos.Offset
                 } deriving (Show)

confMinLength :: Conf -> Int
confMaxLength :: Conf -> Int
confMinLength = fst . confLengths
confMaxLength = snd . confLengths

argConf :: Term Conf
argConf = Conf <$>
          argBamInput <*>
          argOutput <*>
          argBedFiles <*>
          argGeneFiles <*>
          argFlank <*>
          argBody <*>
          argLengths <*>
          argAnnotate <*>
          argBinSize

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

argGeneFiles :: Term [FilePath]
argGeneFiles = value $ optAll [] $ (optInfo ["g", "genes"])
  { optName = "GENES.TXT", optDoc = "Tab-delimited table of Transcript<TAB>Gene (or just Transcript to suppress a transcript)" }

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

argBinSize :: Term Pos.Offset
argBinSize = value $ opt 100000 $ (optInfo ["z", "binsize"])
  { optName = "BINSIZE", optDoc = "Bin size for transcript lookup map" }

main :: IO ()
main = run ( fpframe, info )
  where fpframe = doFpFraming <$> argConf
        info = defTI { termName = "fp-framing"
                     , version = "150304"
                     , termDoc = "Calculates ribosome profiling QC information including reading frame bias and start and stop codon meta-genes"
                     , man = map P [ "Calculates quality control statistics from ribosome profiling data. These QC information are reported in four distinct data files that summarize results from mapping footprint alignments onto protein-coding gene annotations. Individual footprints can also be annotated with QC classifications in a BAM file output."
                                   , "Two files contain a 2-D metagene analysis of footprints around the start and the end of protein-coding genes. These meta-genes include all footprints whose length falls within the range specified by \"lengths\". The \"flanking\" argument indicates the range of positions included in the metagene. For each position (listed in the first column), the second column gives the total number of footprints starting at that position, and the following columns give the number at each individual length, increasing. In the start metagene, written in OUTPUT_around_start.txt, 0 is the first nucleotide of the start codon. In the end metagene, written in OUTPUT_around_end.txt, 0 is the _first_ nucleotide of the stop codon (a change from earlier programs!)."
                                   , "One file contains a table of reading frame position for footprints falling within the CDS of protein-coding genes. The region within the CDS is given by the \"cdsbody\" argument and ensures that the start of the footprint is at least AFTERSTART past the start codon while the end of the footprint is at least BEFOREEND from the stop codon. Each row provides statistics for one footprint length within the range of lengths in the \"lengths\" argument. The \"fract\" column indicates the fraction of all considered footprints at that specific length; the \"N0\", \"N1\", and \"N2\" columns give the count of footprints whose 5\' end is on the first, second, or third nucleotide of a codon; the \""
                                   ]
                     }
