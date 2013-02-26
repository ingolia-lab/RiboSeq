{-# LANGUAGE OverloadedStrings #-}

module Main
       where 

import Control.Applicative
import Control.Exception
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

import qualified Bio.SamTools.Bam as Bam
import qualified Bio.SamTools.BamIndex as BamIndex

import Bio.SeqLoc.Bed
import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.LocRepr
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.Position as Pos
import qualified Bio.SeqLoc.SpliceLocation as SpLoc
import Bio.SeqLoc.Strand
import Bio.SeqLoc.Transcript

import Bio.RiboSeq.BamFile
import Bio.RiboSeq.CodonAssignment
import Bio.RiboSeq.Counting
--import qualified Bio.RiboSeq.QExpr as QE

main :: IO ()
main = getArgs >>= handleOpt . getOpt RequireOrder optDescrs
    where handleOpt (_,    _,         errs@(_:_)) = usage (unlines errs)
          handleOpt (args, [bam], []) = either usage (doFpQuantify bam) $ argsToConf args
          handleOpt (_,    _,     []) = usage "Specify just one sorted, indexed BAM file"
          usage errs = do prog <- getProgName
                          let progline = prog ++ " [OPTIONS] <BAM>"
                          hPutStr stderr $ usageInfo progline optDescrs
                          hPutStrLn stderr errs

doFpQuantify :: FilePath -> Conf -> IO ()
doFpQuantify bam conf = do asite <- readASiteDelta $ confASite conf
                           withFile (confOutput conf) WriteMode $ \hout ->
                             bracket (BamIndex.open bam) BamIndex.close $ \bidx -> 
                             mapOverTranscripts (confBeds conf) $ \trx -> 
                             writeQuant conf asite hout bidx trx

codonLength :: (Loc.Location l) => l -> Int
codonLength l = fromIntegral $ (Loc.length l + 1) `div` 3

writeQuant :: Conf -> ASiteDelta -> Handle -> BamIndex.IdxHandle -> Transcript -> IO ()
writeQuant conf asite h bidx trx = do 
  ct <- newEnumXountIO (False, True)
  let bamLoc = liftM (stranded (confStrand conf)) . Bam.refSpLoc
      countIn = countBam conf ct $! maybe False (isIn conf asite trx) . bamLoc
  _ <- mapOverBams bidx countIn trx
  n <- readEnumXount ct True
  BS.hPutStrLn h $ BS.intercalate "\t" [ unSeqLabel . trxId $ trx
                                       , BS.pack $ featureLengthStr conf trx
                                       , BS.pack $ showFFloat (Just 0) n ""
                                       ]
  hFlush h
  when (confDebug conf) $ do
    readEnumXount ct True >>= print
    readEnumXount ct False >>= print
  where trxname = unSeqLabel . trxId $ trx

featureLengthStr :: Conf -> Transcript -> String
featureLengthStr conf | confWhole conf = trxLengthStr
                      | otherwise      = cdsLengthStr
  where trxLengthStr trx = show . codonLength . unOnSeq . location $! trx
        cdsLengthStr trx = maybe noCdsStr cdsStr $! cds trx
          where noCdsStr = "N/A"
                cdsStr cdsloc = let (inset5, inset3) = confCdsInsets conf
                                    len = codonLength cdsloc - (inset5 + inset3)
                                in show len

isIn :: Conf -> ASiteDelta -> Transcript -> SpLoc.SpliceLoc -> Bool
isIn conf | confWhole conf = isInTrx conf
          | otherwise      = isInCds conf

isInCds :: Conf -> ASiteDelta -> Transcript -> SpLoc.SpliceLoc -> Bool
isInCds conf asite trx = maybe (const False) incds $! cds trx
  where incds cdsloc = let (inset5, inset3) = confCdsInsets conf
                           inbody codon = codon >= inset5 && codon < (codonLength cdsloc - inset3)
                       in \loc -> case cdsReadASite asite trx loc of
                         (ReadASite cdsa) -> inbody . fromIntegral . cfCodon . ntOffsetToCodonFrame $ cdsa
                         _ -> False

isInTrx :: Conf -> ASiteDelta -> Transcript -> SpLoc.SpliceLoc -> Bool
isInTrx conf asite trx loc = case trxReadASite asite trx loc of
  (ReadASite _pos) -> True
  _ -> False

countBam :: (Show a, Enum a) => Conf -> EnumXountIO a -> (Bam.Bam1 -> a) -> Bam.Bam1 -> IO ()
countBam conf count f bam 
  = let x | confNHits conf = maybe 0.0 (recip . fromIntegral) $! Bam.nHits bam
          | otherwise = 1.0
    in xountEnum count (f bam) x

data Conf = Conf { confOutput :: !(FilePath) 
                 , confBeds :: ![FilePath]
                 , confASite :: !FilePath
                 , confCdsInsets :: !(Int, Int)
                 , confStats :: !(Maybe FilePath)
                 , confFraming :: !(Maybe (Int, Int))
                 , confNHits :: !Bool
                 , confWhole :: !Bool
                 , confDebug :: !Bool
                 , confStrand :: !Strand
                 } deriving (Show)

defaultCdsInsets :: (Int, Int)
defaultCdsInsets = (15, 5)

confFramingOutput :: Conf -> FilePath
confFramingOutput conf = (stem ++ "_fr") <.> ext
  where (stem, ext) = case splitExtension $ confOutput conf of
          (stem0, "") -> (stem0, "txt")
          (stem0, ext0) -> (stem0, ext0)

data Arg = ArgOutput { unArgOutput :: !String }
         | ArgBed { unArgBed :: !String }
         | ArgASite { unArgASite :: !String }
         | ArgCdsInsets { unArgCdsInsets :: !String }
         | ArgStats { unArgStats :: !String }
         | ArgFraming { unArgFraming :: !String }
         | ArgNHits
         | ArgWhole
         | ArgDebug
         | ArgReverse
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

argCdsInsets :: Arg -> Maybe String
argCdsInsets (ArgCdsInsets cdsInsets) = Just cdsInsets
argCdsInsets _ = Nothing

argStats :: Arg -> Maybe String
argStats (ArgStats stats) = Just stats
argStats _ = Nothing

argFraming :: Arg -> Maybe String
argFraming (ArgFraming framing) = Just framing
argFraming _ = Nothing

optDescrs :: [OptDescr Arg]
optDescrs = [ Option ['o'] ["output"]  (ReqArg ArgOutput "OUTFILE")      "Output filename"
            , Option ['b'] ["bed"]     (ReqArg ArgBed "BED")             "Bed filename"
            , Option ['a'] ["asite"]   (ReqArg ArgASite "ASITEFILE")     "A site offsets filename"
            , Option ['c'] ["cds-inset"] (ReqArg ArgCdsInsets "IN5',IN3'") "CDS insets for quantitation"
            , Option ['s'] ["stats"]   (ReqArg ArgStats "STATBASE")      "Base filename for statistics"            
            , Option ['f'] ["framing"] (ReqArg ArgFraming "MINLEN,MAXLEN") "Per-gene framing for length range"
            , Option ['n'] ["nhits"]   (NoArg ArgNHits)                    "Scale by NH of alignment"
            , Option ['w'] ["whole"]   (NoArg ArgWhole)                    "Quantify over the whole transcript, not just the CDS"
            , Option ['d'] ["debug"]   (NoArg ArgDebug)                    "Debug read -> feature assignment"
            , Option ['r'] ["reverse"] (NoArg ArgReverse)                "Reverse strand reads"
            ]

argsToConf :: [Arg] -> Either String Conf
argsToConf = runReaderT conf
    where conf = Conf <$> 
                 findOutput <*>
                 findBeds <*>
                 findASite <*>
                 findCdsInsets <*>
                 findStats <*>
                 findFraming <*>
                 ReaderT (return . elem ArgNHits) <*>
                 ReaderT (return . elem ArgWhole) <*>
                 ReaderT (return . elem ArgDebug) <*>
                 findStrand
          findOutput = ReaderT $ maybe (Left "No out base") return  . listToMaybe . mapMaybe argOutput
          findBeds = ReaderT $ return . mapMaybe argBed
          findASite = ReaderT $ maybe (Left "No A sites") return . listToMaybe . mapMaybe argASite
          findCdsInsets = ReaderT $ maybe (return defaultCdsInsets) parseIntPair . listToMaybe . mapMaybe argCdsInsets
          findStats = ReaderT $ return . listToMaybe . mapMaybe argStats
          findFraming = ReaderT $ maybe (return Nothing) (liftM Just . parseIntPair) . listToMaybe . mapMaybe argFraming
          findStrand = ReaderT $ \args -> if elem ArgReverse args then return Minus else return Plus
                        
parseIntPair :: (Integral a) => String -> Either String (a, a)
parseIntPair = AP.parseOnly intPair . BS.pack
  where intPair = (,) <$> AP.signed AP.decimal <*> 
                  (AP.char ',' *> AP.signed AP.decimal <* AP.endOfInput)
                  
debugPut :: Conf -> String -> IO ()
debugPut conf | confDebug conf = putStrLn
              | otherwise      = const $ return ()
                                 