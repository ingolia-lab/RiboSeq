module Bio.RiboSeq.StartSite
       where

import Control.Arrow
import Control.Monad
import qualified Data.ByteString.Char8 as BS
import Data.List (intercalate)
import Data.Monoid
import Numeric
import System.Directory
import System.Exit
import System.FilePath

import qualified Bio.SamTools.BamIndex as BamIndex
import qualified Data.Conduit as C
import qualified Data.Conduit.List as C
import qualified Data.Vector.Unboxed as U

import Bio.SeqLoc.Transcript

import Bio.RiboSeq.BamFile
import Bio.RiboSeq.CodonAssignment
import Bio.RiboSeq.SVMLight
import Bio.RiboSeq.StartSVM

data LtmSamples  = LtmSamples  { ltmSample, chxSample :: !HarrSample } deriving (Read, Show)
data LtmModel = LtmModel { ltmMinReads :: !Int
                         , ltmMinDiff :: !Double
                         , ltmMinRatio :: !Double
                         } deriving (Read, Show)

defaultLtmModel :: LtmModel
defaultLtmModel = LtmModel { ltmMinReads = 10
                           , ltmMinDiff = 5.0
                           , ltmMinRatio = 2.0
                           }

data StartModel = StartModel { harrModel :: !HarrModel
                             , ltmSamples :: !LtmSamples
                             , ltmModel :: !LtmModel
                             } deriving (Read, Show)

data TrxLtmProf = TrxLtmProf { ltmProf, chxProf :: !TrxProfile }

chxTotal :: TrxLtmProf -> Int
chxTotal = U.sum . trxVectorV . chxProf

readLtmProfile :: Transcript -> LtmSamples -> IO TrxLtmProf
readLtmProfile trx (LtmSamples ltmsample chxsample) = do 
  ltmprof <- readHarrProfile trx ltmsample
  chxprof <- readHarrProfile trx chxsample
  return $! TrxLtmProf ltmprof chxprof

testLtm :: LtmModel -> LtmSamples -> TrainPosns -> [Transcript] -> IO TestScore
testLtm model samples posns trxset = trxsrc C.$$ C.foldMap testLtmTrx
  where trxsrc = C.sourceList trxset C.$= C.mapM (flip readLtmProfile samples)
        testLtmTrx prof = maybe mempty testLtmCds $ trxVectorStart . ltmProf $ prof
          where testLtmCds cdsstart
                  | chxTotal prof < 1 = mempty
                  | otherwise = let score = ((ltmIsStart model prof &&& ltmStartRatio prof) . (+ cdsstart)) 
                                    starts = map score (posnStarts posns)
                                    nonstarts = map score (posnNonStarts posns)
                        in TestData { truePos = U.fromList . map snd . filter fst $ starts
                                    , falsePos = U.fromList . map snd . filter fst $ nonstarts
                                    , falseNeg = U.fromList . map snd . filter (not . fst) $ starts
                                    , trueNeg = U.fromList . map snd . filter (not . fst) $ nonstarts
                                    }
        
ltmStartProfile :: LtmModel -> TrxLtmProf -> TrxVector Bool
ltmStartProfile model prof = TrxVector (trxVectorTrx . ltmProf $ prof) stvec
  where stvec = U.generate (trxVectorLength . ltmProf $ prof) (ltmIsStart model prof)
        
ltmIsStart :: LtmModel -> TrxLtmProf -> Int -> Bool
ltmIsStart model (TrxLtmProf (TrxVector trx ltm) (TrxVector _ chx)) st
  | st + 4 >= U.length ltm = False
  | otherwise = let ltmc = fromIntegral . U.sum . U.take 3 . U.drop (st + 2) $ ltm
                    chxc = fromIntegral . U.sum . U.take 3 . U.drop (st + 2) $ ltm
                    ltmttl = fromIntegral . U.sum $ ltm
                    chxttl = fromIntegral . U.sum $ chx
                    expc = chxc * ltmttl / chxttl
                in and [ chxttl > 0.5
                       , ltmc >= fromIntegral (ltmMinReads model)
                       , ltmc - expc >= ltmMinDiff model
                       , ltmc >= expc * ltmMinRatio model
                       ]

ltmStartRatio :: TrxLtmProf -> Int -> Double
ltmStartRatio (TrxLtmProf (TrxVector trx ltm) (TrxVector _ chx)) st
  | st + 4 >= U.length ltm = 0.0
  | otherwise = let ltmc = fromIntegral . U.sum . U.take 3 . U.drop (st + 2) $ ltm
                    chxc = fromIntegral . U.sum . U.take 3 . U.drop (st + 2) $ ltm
                    ltmttl = fromIntegral . U.sum $ ltm
                    chxttl = fromIntegral . U.sum $ chx
                    expc = chxc * ltmttl / chxttl
                in if chxttl < 0.5 then -1.0 else (if chxc < 0.5 then ltmc else ltmc / expc)

ltmTestScore :: TrxLtmProf -> Int -> String
ltmTestScore (TrxLtmProf ltmp@(TrxVector _ ltm) chxp@(TrxVector _ chx)) dstart
  = intercalate "\t" . addname . maybe ["n/a"] (atstart . (+ dstart)) $ trxVectorStart ltmp
  where addname = ((BS.unpack $ trxVectorName ltmp) :)
        atstart st 
          | st < 0 = [ "n/a" ]
          | otherwise = let ltm_st = ltm U.! (st + 3)
                            ltm_ttl = fromIntegral $ U.sum ltm
                            ltm_max | (st - 1 >= 0) = (U.maxIndex . U.take 9 . U.drop (st - 1) $ ltm) - 1
                                    | otherwise     = (U.maxIndex . U.take 7                   $ ltm)
                            chx_st = chx U.! (st + 3)
                            chx_ttl = fromIntegral $ U.sum chx
                            showz :: Double -> String
                            showz x = showFFloat (Just 0) x ""
                            showf :: Double -> String
                            showf x = showFFloat (Just 4) x ""
                            ltm_c = U.sum . U.take 3 . U.drop (st + 2) $ ltm
                            chx_c = U.sum . U.take 3 . U.drop (st + 2) $ chx
                            exp_c | chx_ttl > 0.5 = (fromIntegral chx_c) * ltm_ttl / chx_ttl
                                  | otherwise = 0.0
                        in [ show ltm_st, show chx_st, showz ltm_ttl, showz chx_ttl
                           , showf (fromIntegral ltm_st / ltm_ttl)
                           , showf (fromIntegral chx_st / chx_ttl)
                           , show ltm_max
                           , show ltm_c, show chx_c
                           , showz exp_c, showz (fromIntegral ltm_c - exp_c), showf (fromIntegral ltm_c / exp_c)
                           ]

-- Identiﬁcation of TIS Positions. A peak is deﬁned at the nucleotide
-- level on a transcript. A peak position satisﬁes the following conditions: (i) The transcript must have both LTM and CHX reads.
-- (ii) The position must have at least 10 reads from the LTM data.
-- (iii) The position must be a local maximum within seven nucleotides (4). The position must have RLTM − RCHX of at least 0.05,
-- where Rk = (Xk/Nk) × 10 (k = LTM, CHX), Xk is the number of
-- reads on that position in data k, and Nk
-- is the total number of reads
-- on that transcript in data k. Generally, a peak position is also
-- designated a TIS. However, if a peak was not detected on the ﬁrst
-- position of any AUG or near-cognate start codon but was present
-- at the ﬁrst position of a codon immediately preceding or succeeding one of these codons, the position was designated a TIS.

-- Assume 1000 nt CDS
-- Difference (LTM - EXP) = (Xltm - (Xchx * Nltm / Nchx)) = 0.005 * Nltm
-- Ratio (LTM / EXP) 
--   Rltm - Rchx >= 0.05
--   Rltm >= 0.05 + Rchx
--   