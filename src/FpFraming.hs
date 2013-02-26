{-# LANGUAGE OverloadedStrings #-}

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
import System.Console.GetOpt
import System.Environment
import System.IO

import qualified Data.Attoparsec as AP
import qualified Data.Attoparsec.Char8 as AP
import qualified Data.Vector as V

import qualified Bio.SamTools.BamIndex as BamIndex
import Bio.SeqLoc.LocRepr
import qualified Bio.SeqLoc.Position as Pos

import Bio.RiboSeq.BamFile
import Bio.RiboSeq.Framing

main :: IO ()
main = getArgs >>= handleOpt . getOpt RequireOrder optDescrs
    where handleOpt (_,    _,         errs@(_:_)) = usage (unlines errs)
          handleOpt (args, bams, []) = either usage (doFpFraming bams) $ argsToConf args
          usage errs = do prog <- getProgName
                          let progline = prog ++ " [OPTIONS] <BAM1> <BAM2> ..."
                          hPutStr stderr $ usageInfo progline optDescrs
                          hPutStrLn stderr errs

doFpFraming :: [FilePath] -> Conf -> IO ()
doFpFraming bams conf = do trstartio <- newTerminus (confFlank conf) (confLengths conf)
                           trendio <- newTerminus (confFlank conf) (confLengths conf)
                           frameio <- newFrame (confLengths conf)
                           let count trx iloc = do countAtStart trstartio trx iloc
                                                   countAtEnd trendio trx iloc
                                                   countInCds (confCdsBody conf) frameio trx iloc
                           withMany (\bam -> bracket (BamIndex.open bam) (BamIndex.close)) bams $ \bidxs ->
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

data Conf = Conf { confOutput :: !(FilePath) 
                 , confBeds :: [FilePath]
                 , confFlank :: !(Pos.Offset, Pos.Offset)
                 , confCdsBody :: !(Pos.Offset, Pos.Offset)
                 , confLengths :: !(Int, Int)
                 } deriving (Show)

defaultFlank :: (Pos.Offset, Pos.Offset)
defaultFlank = (-100, 100)

defaultCdsBody :: (Pos.Offset, Pos.Offset)
defaultCdsBody = (34, 31)

defaultLengths :: (Int, Int)
defaultLengths = (25, 34)

data Arg = ArgOutput { unArgOutput :: !String }
         | ArgBed { unArgBed :: !String }
         | ArgFlanking { unArgFlanking :: !String }
         | ArgCdsBody { unArgCdsBody :: !String }
         | ArgLengths { unArgLengths :: !String }
         deriving (Show, Read, Eq, Ord)

argOutput :: Arg -> Maybe String
argOutput (ArgOutput del) = Just del
argOutput _ = Nothing

argBed :: Arg -> Maybe String
argBed (ArgBed bed) = Just bed
argBed _ = Nothing

argFlanking :: Arg -> Maybe String
argFlanking (ArgFlanking flanking) = Just flanking
argFlanking _ = Nothing

argCdsBody :: Arg -> Maybe String
argCdsBody (ArgCdsBody cdsbody) = Just cdsbody
argCdsBody _ = Nothing

argLengths :: Arg -> Maybe String
argLengths (ArgLengths lengths) = Just lengths
argLengths _ = Nothing

optDescrs :: [OptDescr Arg]
optDescrs = [ Option ['o'] ["output"]   (ReqArg ArgOutput "OUTFILE")        "Output filename"
            , Option ['b'] ["bed"]      (ReqArg ArgBed "BED")               "Bed filename"
            , Option ['f'] ["flanking"] (ReqArg ArgFlanking "START,END")    flankingDesc
            , Option ['c'] ["cdsbody"]  (ReqArg ArgCdsBody "5',3'")         cdsbodyDesc
            , Option ['l'] ["lengths"]  (ReqArg ArgLengths "MIN,MAX")       lengthDesc
            ]
  where flankingDesc = "Range of profile before and after terminus [" ++
                       (BS.unpack . repr . fst $ defaultFlank) ++ 
                       "," ++ (BS.unpack . repr . snd $ defaultFlank) ++ "]"
        cdsbodyDesc = "Inset for footprint positions in the CDS body [" ++
                      (BS.unpack . repr . fst $ defaultCdsBody) ++ 
                      "," ++ (BS.unpack . repr . snd $ defaultCdsBody) ++ "]"
        lengthDesc = "Footprint fragment length range [" ++
                     (show . fst $ defaultLengths) ++ 
                     "," ++ (show . snd $ defaultLengths) ++ "]"
        
argsToConf :: [Arg] -> Either String Conf
argsToConf = runReaderT conf
    where conf = Conf <$> 
                 findOutput <*>
                 findBeds <*>
                 findFlank <*>
                 findCdsBody <*>
                 findLengths
          findOutput = ReaderT $ maybe (Left "No out base") return  . listToMaybe . mapMaybe argOutput
          findBeds = ReaderT $ return . mapMaybe argBed
          findFlank = ReaderT $ maybe (return defaultFlank) parseIntPair . listToMaybe . mapMaybe argFlanking
          findCdsBody = ReaderT $ maybe (return defaultCdsBody) parseIntPair . listToMaybe . mapMaybe argCdsBody
          findLengths = ReaderT $ maybe (return defaultLengths) parseIntPair . listToMaybe . mapMaybe argLengths
          
parseIntPair :: (Integral a) => String -> Either String (a, a)
parseIntPair = AP.parseOnly intPair . BS.pack
  where intPair = (,) <$> AP.signed AP.decimal <*> 
                  (AP.char ',' *> AP.signed AP.decimal <* AP.endOfInput)