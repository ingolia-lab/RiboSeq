{-# LANGUAGE BangPatterns, OverloadedStrings #-}
module Main
       where

import Control.Applicative
import Control.Monad
import Control.Monad.Reader
import Control.Monad.Error
import qualified Data.ByteString.Char8 as BS
import qualified Data.HashSet as HS
import Data.List
import Data.Maybe
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U

import Numeric
import System.Console.CmdTheLine
import System.Environment
import System.IO

import qualified Data.Attoparsec.ByteString.Char8 as AP

import Statistics.Test.MannWhitneyU
import Statistics.Quantile

main :: IO ()
main = run ( subset, info )
  where info = defTI { termName = "subset-mw"
                     , version = "0.0"
                     , termDoc = "Mann-Whitney U test for GO enrichment in ordinal data"
                     }
        subset = subsetMannWhitney <$> argConf

subsetMannWhitney :: Conf -> IO ()
subsetMannWhitney conf = do
  dvals <- readData conf
  subsets <- readSubsets conf
  let stats = V.fromList $ mapMaybe (subsetStat conf dvals) subsets
  writeFile (cOutBase conf ++ "_stats.txt") $ unlines . map show . V.toList $ stats
  let (fdrEsts, fdrTrials) = statFdr conf stats
  BS.writeFile (cOutBase conf ++ "_fdr_est.txt") $ fdrTrialReport conf fdrTrials
  BS.writeFile (cOutBase conf ++ "_subsets.txt") $ subsetReport conf stats fdrEsts
           
subsetReport :: Conf -> V.Vector Stat -> V.Vector FdrEst -> BS.ByteString
subsetReport _conf stats fdrs = BS.unlines $ header : subsetLines
  where header = BS.intercalate "\t" [ "# GO", "MedDiff", "MedIn", "MedOut", "p", "q", "NumIn" ]
        subsetLines = V.toList $ V.zipWith subsetLine stats fdrs
        subsetLine stat fdr = BS.pack . intercalate "\t" $
                              [ BS.unpack . sName . sSubset $ stat
                              , showFFloat (Just 2) (sMedDiff stat) ""
                              , showFFloat (Just 2) (sMedIn stat) ""
                              , showFFloat (Just 2) (sMedOut stat) ""
                              , showEFloat (Just 1) (estP fdr) ""
                              , showEFloat (Just 1) (estQ fdr) ""
                              , show . HS.size . sKeys . sSubset $ stat
                              ]

fdrTrialReport :: Conf -> [FdrTrial] -> BS.ByteString
fdrTrialReport _conf trials = BS.unlines $ header : map trialLine trials
  where header = BS.intercalate "\t" [ "# p", "N_tot", "N_sig", "Exp_false", "Qest" ]
        trialLine trial = BS.pack . intercalate "\t" $
                          [ showEFloat (Just 1) (fdrP trial) ""
                          , show $ fdrNttl trial
                          , show $ fdrNsig trial
                          , showFFloat (Just 1) ((fromIntegral $ fdrNttl trial) * (fdrP trial)) ""
                          , showEFloat (Just 1) (fdrQ trial) ""
                          ]

data FdrTrial = FdrTrial { fdrP :: !Double
                         , fdrQ :: !Double
                         , fdrNttl :: !Int
                         , fdrNsig :: !Int
                         , fdrSigs :: !(U.Vector Bool)
                         } deriving (Show)

mannWhitneyFdrTrial :: Double -> V.Vector Stat -> FdrTrial
mannWhitneyFdrTrial p stats = FdrTrial p q nttl nsig sigs
  where issig s = maybe False (== Significant) $ mannWhitneyUSignificant TwoTailed (sNIn s, sNOut s) p (sUIn s, sUOut s)
        sigs = U.convert . V.map issig $ stats
        nttl = V.length stats
        nsig = U.foldl' (\n s -> if s then succ n else n) 0 $ sigs
        q = p * (fromIntegral nttl) / (fromIntegral nsig)

fdrTable :: Conf -> V.Vector Stat -> [FdrTrial]
fdrTable conf stats = take maxPCount $ unfoldr genFdr $ cPMax conf
  where genFdr p | p < (cPMin conf) = Nothing
                 | otherwise = case mannWhitneyFdrTrial p stats of
                   trial | U.or (fdrSigs trial) -> Just (trial, p / cPStep conf)
                         | otherwise -> Nothing
        maxPCount = 50

data FdrEst = FdrEst { estP, estQ :: !Double } deriving (Show)

statFdr :: Conf -> V.Vector Stat -> (V.Vector FdrEst, [FdrTrial])
statFdr conf stats = let !fdrests = V.imap estfdr $ stats
                     in (fdrests, trials)
  where trials = fdrTable conf stats
        estfdr i _s = let sigtrials = filter (\trial -> (fdrSigs trial) U.! i)  trials
                          sigps = map fdrP sigtrials
                          sigqs = map fdrQ sigtrials
                      in FdrEst (minimum $ 1.0 : sigps) (minimum $ 1.0 : sigqs)

data Stat = Stat { sSubset :: !Subset
                 , sNIn, sNOut :: !Int
                 , sUIn, sUOut :: !Double
                 , sMedIn, sMedOut :: !Double
                 } deriving (Show)

sMedDiff :: Stat -> Double
sMedDiff s = sMedIn s - sMedOut s

subsetStat :: Conf -> [Datum] -> Subset -> Maybe Stat
subsetStat conf dats subset
  | min ni no >= fromMaybe 1 (cMinSize conf) = Just $! Stat subset ni no ui uo medi medo
  | otherwise = Nothing
  where (vin, vout) = partitionData (flip HS.member $ sKeys subset) dats
        ni = U.length vin
        no = U.length vout
        medi = continuousBy medianUnbiased 1 2 vin
        medo = continuousBy medianUnbiased 1 2 vout
        (ui, uo) = mannWhitneyU vin vout

partitionData :: (BS.ByteString -> Bool) -> [Datum] -> (U.Vector Double, U.Vector Double)
partitionData isin dats = (tovec dvin, tovec dvout)
  where (dvin, dvout) = partition (isin . dKey) dats
        tovec = U.fromList . map dValue

data Subset = Subset { sName :: !BS.ByteString
                     , sKeys :: !(HS.HashSet BS.ByteString)
                     } deriving (Show)

readSubsets :: Conf -> IO [Subset]
readSubsets conf = (BS.readFile . cSubsetFile $ conf) >>= parseSubsets
  where parseSubsets = mapM parseSubset . BS.lines
        parseSubset l = case BS.split '\t' l of
          (key:vals@(_:_)) -> return $! Subset key (HS.fromList vals)
          _ -> ioError . userError $ unlines [ "Error reading " ++ show (cSubsetFile conf) ++ ":\n  Malformed line " ++ show l ]

data Datum = Datum { dKey :: !BS.ByteString
                   , dValue :: !Double
                   } deriving (Show)

readData :: Conf -> IO [Datum]
readData conf = (BS.readFile . cDataFile $ conf) >>= parseData
  where parseData = mapM parseDatum . BS.lines
        parseDatum l = case BS.split '\t' l of
          [key, valstr] -> case AP.parseOnly AP.double valstr of
            (Right val) -> return $! Datum key val
            (Left _e) -> err $ "Malformed value in " ++ show l
          _ -> err $ "Malformed data line " ++ show l         
        err errstr = ioError . userError $ unlines [ "Reading " ++ show (cDataFile conf) ++ ":\n  " ++ errstr ] 

data Conf = Conf { cSubsetFile :: !FilePath
                 , cDataFile :: !FilePath
                 , cOutBase :: !FilePath
                 , cMinSize :: !(Maybe Int)
                 , cPMin :: !Double
                 , cPMax :: !Double
                 , cPStep :: !Double
                 } deriving (Show)

argConf :: Term Conf
argConf = Conf <$>
          argSubsetFile <*>
          argDataFile <*>
          argOutBase <*>
          argMinSize <*>
          pure 1e-9 <*>
          pure 1e-1 <*>
          pure 2.0

argSubsetFile :: Term FilePath
argSubsetFile = required $ opt Nothing $ (optInfo ["s", "subsets"])
  { optName = "SUBSETS.TXT", optDoc = "Table of subsets with lines of subset<TAB>gene1,gene2,...,geneN" }

argDataFile :: Term FilePath
argDataFile = required $ opt Nothing $ (optInfo ["d", "data"])
  { optName = "DATA.TXT", optDoc = "Table of data with lines of gene<TAB>datum" }

argOutBase :: Term FilePath
argOutBase = required $ opt Nothing $ (optInfo ["o", "outbase"])
  { optName = "OUTBASE", optDoc = "Base name for constructing output filenames" }

argMinSize :: Term (Maybe Int)
argMinSize = value $ opt Nothing $ (optInfo ["n", "min-size"])
  { optName = "MIN-SIZE", optDoc = "Minimum number of genes inside and outside a subset for comparison" }
