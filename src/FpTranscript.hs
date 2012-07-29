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

import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U

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

main :: IO ()
main = getArgs >>= handleOpt . getOpt RequireOrder optDescrs
    where handleOpt (_,    _,         errs@(_:_)) = usage (unlines errs)
          handleOpt (args, [bam], []) = either usage (doFpTranscript bam) $ argsToConf args
          handleOpt (_,    _,     []) = usage "Specify just one sorted, indexed BAM file"
          usage errs = do prog <- getProgName
                          let progline = prog ++ " [OPTIONS] <BAM>"
                          hPutStr stderr $ usageInfo progline optDescrs
                          hPutStrLn stderr errs

doFpTranscript :: FilePath -> Conf -> IO ()
doFpTranscript bam conf = do asites <- readASite $ confASite conf
                             bracket (BamIndex.open bam) BamIndex.close $ \bidx ->
                               findTranscript (confBeds conf) (confTranscript conf) >>= \trx ->
                               doTranscript conf asites bidx trx

doTranscript :: Conf -> ASites -> BamIndex.IdxHandle -> Transcript -> IO ()
doTranscript conf
  | confByLength conf  = doTranscriptByLength conf
  | confCdsCodons conf = doTranscriptCodon conf
  | otherwise          = doTranscriptNt conf
                        
doTranscriptByLength :: Conf -> ASites -> BamIndex.IdxHandle -> Transcript -> IO ()
doTranscriptByLength conf asites bidx trx
  = do prof <- transcriptNtLengthProfile asites bidx trx
       msequ <- maybe (return Nothing) (flip FaIdx.readLoc (location trx)) (confFasta conf)
       let (hdrs, pf) = ntLengthProfile trx msequ prof
       writeFile (confOutput conf) pf
       writeInfo conf hdrs

doTranscriptCodon :: Conf -> ASites -> BamIndex.IdxHandle -> Transcript -> IO ()
doTranscriptCodon conf asites bidx trx
  = do (hdrs, pf) <- liftM (codonProfile trx) $ transcriptNtProfile (aSiteDelta asites) bidx trx
       writeInfo conf hdrs
       writeFile (confOutput conf) pf

doTranscriptNt :: Conf -> ASites -> BamIndex.IdxHandle -> Transcript -> IO ()
doTranscriptNt conf asites bidx trx
  = do prof <- transcriptNtProfile (aSiteDelta asites) bidx trx
       msequ <- maybe (return Nothing) (flip FaIdx.readLoc (location trx)) (confFasta conf)
       let (hdrs, pf) = ntProfile trx msequ prof
       writeFile (confOutput conf) pf
       writeInfo conf hdrs

ntProfile :: Transcript -> Maybe BS.ByteString -> U.Vector Int -> (String, String)
ntProfile trx msequ prof = (unlines headers, unlines body)
  where headers = map ("# " ++ ) [ BS.unpack . unSeqLabel . trxId $ trx
                                 , reprStr . location $ trx
                                 , junctionLocs trx
                                 , maybe "n/a" reprStr . cds $ trx                                 
                                 ] ++ maybe [] (flip ntStats prof) (cds trx)
        body = map profline [0..(U.length prof - 1)]
        profline i = intercalate "\t" fields
          where fields = profileLocFields trx i ++ maybeNt ++ maybeAa ++ [ totalfield ]
                totalfield = show $ prof U.! i
                maybeNt = maybe [] ntAt msequ
                maybeAa = maybe [] aaAt msequ
                  where aaAt sequ = [ indent $ maybe "x" show $ lookup (BS.take 3 $ BS.drop i sequ) trans_tbl ]
                        indent = ((replicate (2 * (i `mod` 3)) ' ') ++)
                ntAt sequ = [ [BS.index sequ i] ]

codonProfile :: Transcript -> U.Vector Int -> (String, String)
codonProfile trx prof = (unlines headers, unlines $ maybe [] body $ cds trx)
  where headers = map ("# " ++ ) [ BS.unpack . unSeqLabel . trxId $ trx
                                 , reprStr . location $ trx
                                 , maybe "n/a" reprStr . cds $ trx
                                 ] ++ maybe [] (flip ntStats prof) (cds trx)
        body cdsloc = map profline [0..(U.length cdsprof - 1)]
          where cdsprof = profileFromStart Nothing cdsloc prof
                profline i = intercalate "\t" [ show i, show (cdsprof U.! i) ]

ntLengthProfile :: Transcript -> Maybe BS.ByteString -> V.Vector (EnumCount Pos.Offset) -> (String, String)
ntLengthProfile trx msequ prof = (unlines headers, unlines body)
  where headers = map ("# " ++ ) $ [ BS.unpack . unSeqLabel . trxId $ trx
                                   , reprStr . location $ trx
                                   , junctionLocs trx
                                   , maybe "n/a" reprStr . cds $ trx
                                   ] ++ maybe [] (flip lengthStats prof) (cds trx)
        body = map profline [0..(V.length prof - 1)]
        profline i = intercalate "\t" fields
          where fields = profileLocFields trx i ++ maybeNt ++ [ totalfield ] ++ lenfields
                ec@(EnumCount lb ub _ct) = prof V.! i
                totalfield = show $ enumCountTotal ec
                lenfields = [ show (indexEnumCount ec l) | l <- [lb..ub] ]
                maybeNt = maybe [] ntAt msequ 
                ntAt sequ = [ [BS.index sequ i] ]

profileLocFields :: Transcript -> Int -> [String]
profileLocFields trx i = [ show i ] ++ cdsFields
  where cdsFields = maybe noCds yesCds $ cds trx
        noCds = [ "n/a", "n/a", "n/a" ]
        yesCds cdsloc = [ show cdsi
                        , show . Pos.unOff . cfCodon $ cfr
                        , show . Pos.unOff . cfFrame $ cfr
                        ]
          where cdsi = i - (fromIntegral . fst . Loc.bounds $ cdsloc)
                cfr = ntOffsetToCodonFrame . fromIntegral $ cdsi

junctionLocs :: Transcript -> String
junctionLocs trx = unwords . map junctionLoc . junctions $ sploc
  where (OnSeq _name sploc) = location trx
        locoff = maybe "***" (reprStr . Pos.offset) . flip Loc.posInto sploc
        junctionLoc j = locoff (donor j) ++ "/" ++ locoff (acceptor j)

lengthStats :: (Loc.Location l) => l -> V.Vector (EnumCount Pos.Offset) -> [String]
lengthStats cdsloc lenprof = [ intercalate "\t" fields ] ++ ntStats cdsloc ttlprof
  where fields = [ "Lengths", reprStr lb ++ "-" ++ reprStr ub ] ++
                 [ show $ totalAt l | l <- [lb..ub] ] ++
                 [ showFFloat (Just 2) (fractAt l) "" | l <- [lb..ub] ]
        (EnumCount lb ub _ct) = lenprof V.! 0
        totalAt l = V.sum . V.map (flip indexEnumCount l) . V.ifilter inCds $ lenprof
        totalAll = (fromIntegral $ sum [ totalAt l | l <- [lb..ub] ]) :: Double
        fractAt l = fromIntegral (totalAt l) / totalAll
        inCds off _ct = maybe False (const True) . ntCodonFrameFromStart cdsloc $ pos
          where pos = (Pos.Pos (fromIntegral off) Plus)
        ttlprof = V.convert . V.map enumCountTotal $ lenprof

ntStats :: (Loc.Location l) => l -> U.Vector Int -> [String]
ntStats cdsloc prof = [ densLine, "Framing\t" ++ profileFraming cdsloc prof ]
  where densLine = intercalate "\t" [ "Density", show ttl, reprStr codonlen, showFFloat (Just 4) dens "" ]
        ttl = U.sum . U.ifilter inCds $ prof
        inCds off _ct = maybe False (const True) . ntCodonFrameFromStart cdsloc $ pos
          where pos = (Pos.Pos (fromIntegral off) Plus)
        codonlen = (Loc.length cdsloc) `div` 3
        dens = (fromIntegral ttl) / (fromIntegral codonlen) :: Double

profileFraming :: (Loc.Location l) => l -> U.Vector Int -> String
profileFraming l prof = intercalate "\t" [ show frneg1, show frzero, show frpos1
                                         , showFFloat (Just 2) fractneg1 ""
                                         , showFFloat (Just 2) fractzero ""
                                         , showFFloat (Just 2) fractpos1 ""
                                         ]
  where frneg1 = U.sum . U.ifilter (hasFrame (-1)) $ prof
        frzero = U.sum . U.ifilter (hasFrame 0) $ prof
        frpos1 = U.sum . U.ifilter (hasFrame 1) $ prof
        total = (fromIntegral $ frneg1 + frzero + frpos1) :: Double
        fractneg1 = fromIntegral frneg1 / total
        fractzero = fromIntegral frzero / total
        fractpos1 = fromIntegral frpos1 / total
        hasFrame fr off _ct = maybe False ((== fr) . cfFrame) . ntCodonFrameFromStart l $ pos
          where pos = (Pos.Pos (fromIntegral off) Plus)

writeInfo :: Conf -> String -> IO ()
writeInfo conf = maybe (const $ return ()) writeFile $ confInfoFile conf

data Conf = Conf { confOutput :: !(FilePath) 
                 , confBeds :: ![FilePath]
                 , confASite :: !FilePath
                 , confTranscript :: !BS.ByteString
                 , confFasta :: !(Maybe FilePath)
                 , confInfo :: !(Maybe (Maybe String))
                 , confCdsCodons :: !Bool
                 , confByLength :: !Bool
                 } deriving (Show)

confInfoFile :: Conf -> Maybe String
confInfoFile conf = liftM (fromMaybe defaultInfo) $ confInfo conf
  where defaultInfo = let (base, ext) = splitExtension $ confOutput conf
                      in addExtension (base ++ "_info") ext

data Arg = ArgOutput { unArgOutput :: !String }
         | ArgBed { unArgBed :: !String }
         | ArgASite { unArgASite :: !String }
         | ArgTranscript { unArgTranscript :: !String }
         | ArgFasta { unArgFasta :: !String }
         | ArgInfo { unArgInfo :: !(Maybe String) }
         | ArgCdsCodons
         | ArgByLength
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

argTranscript :: Arg -> Maybe String
argTranscript (ArgTranscript trx) = Just trx
argTranscript _ = Nothing

argFasta :: Arg -> Maybe String
argFasta (ArgFasta fa) = Just fa
argFasta _ = Nothing

argInfo :: Arg -> Maybe (Maybe String)
argInfo (ArgInfo info) = Just info
argInfo _ = Nothing

optDescrs :: [OptDescr Arg]
optDescrs = [ Option ['o'] ["output"]     (ReqArg ArgOutput "OUTFILE")   "Output filename"
            , Option ['b'] ["bed"]        (ReqArg ArgBed "BED")          "Bed filename"
            , Option ['a'] ["asite"]      (ReqArg ArgASite "ASITEFILE")  "A site offsets filename"
            , Option ['t'] ["transcript"] (ReqArg ArgTranscript "TRX")   "Transcript ID"
            , Option ['c'] ["cds-codons"] (NoArg ArgCdsCodons)           "Codon-indexed CDS profile"
            , Option ['l'] ["lengths"]    (NoArg ArgByLength)            "Length stratified CDS profile"
            , Option ['f'] ["fasta"]      (ReqArg ArgFasta "FASTA")      "Indexed fasta file for sequence"
            , Option ['i'] ["info"]       (OptArg ArgInfo "FILENAME")    "CDS statistics to file"
            ]

argsToConf :: [Arg] -> Either String Conf
argsToConf = runReaderT conf
    where conf = Conf <$> 
                 findOutput <*>
                 findBeds <*>
                 findASite <*>
                 findTranscriptArg <*>
                 findFastaArg <*>
                 findInfoArg <*>
                 ReaderT (return . elem ArgCdsCodons) <*>
                 ReaderT (return . elem ArgByLength)
          findOutput = ReaderT $ maybe (Left "No out base") return  . listToMaybe . mapMaybe argOutput
          findBeds = ReaderT $ return . mapMaybe argBed
          findASite = ReaderT $ maybe (Left "No A sites") return . listToMaybe . mapMaybe argASite
          findTranscriptArg = ReaderT $ maybe (Left "No transcript") (return . BS.pack) . listToMaybe . mapMaybe argTranscript
          findFastaArg = ReaderT $ return . listToMaybe . mapMaybe argFasta
          findInfoArg = ReaderT $ return . listToMaybe . mapMaybe argInfo

data Amino = Ala | Arg | Asn | Asp | Cys | Gln | Glu | Gly
           | His | Ile | Leu | Lys | Met | Phe | Pro | Ser
           | Thr | Tyr | Trp | Val 
           | STP
     deriving (Show,Eq)
              
trans_tbl :: [(BS.ByteString,Amino)]
trans_tbl = [("GCT",Ala),("GCC",Ala),("GCA",Ala),("GCG",Ala),
             ("CGT",Arg),("CGC",Arg),("CGA",Arg),("CGG",Arg),("AGA",Arg),("AGG",Arg),
             ("AAT",Asn),("AAC",Asn),
             ("GAT",Asp),("GAC",Asp),
             ("TGT",Cys),("TGC",Cys),
             ("CAA",Gln),("CAG",Gln),
             ("GAA",Glu),("GAG",Glu),
             ("GGT",Gly),("GGC",Gly),("GGA",Gly),("GGG",Gly),
             ("CAT",His),("CAC",His),
             ("ATT",Ile),("ATC",Ile),("ATA",Ile),
             ("TTA",Leu),("TTG",Leu),("CTT",Leu),("CTC",Leu),("CTA",Leu),("CTG",Leu),
             ("AAA",Lys),("AAG",Lys),
             ("ATG",Met),
             ("TTT",Phe),("TTC",Phe),
             ("CCT",Pro),("CCC",Pro),("CCA",Pro),("CCG",Pro),
             ("TCT",Ser),("TCC",Ser),("TCA",Ser),("TCG",Ser),("AGT",Ser),("AGC",Ser),
             ("ACT",Thr),("ACC",Thr),("ACA",Thr),("ACG",Thr),
             ("TAT",Tyr),("TAC",Tyr),
             ("TGG",Trp),
             ("GTT",Val),("GTC",Val),("GTA",Val),("GTG",Val),
             ("TAG",STP),("TGA",STP),("TAA",STP)
            ]