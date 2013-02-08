module Bio.RiboSeq.StartSite
       where

import qualified Data.ByteString.Char8 as BS
import Data.List (intercalate)
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
data LtmModel = LtmModel { ltmSamples :: !LtmSamples
                         , ltmMinReads :: !Double
                         } deriving (Read, Show)

data StartModel = StartModel { harrModel :: !HarrModel
                             , ltmModel :: !LtmModel
                             } deriving (Read, Show)

data TrxLtmProf = TrxLtmProf { ltmProf, chxProf :: !TrxProfile }

readLtmProfile :: Transcript -> LtmSamples -> IO TrxLtmProf
readLtmProfile trx (LtmSamples ltmsample chxsample) = do 
  ltmprof <- readHarrProfile trx ltmsample
  chxprof <- readHarrProfile trx chxsample
  return $! TrxLtmProf ltmprof chxprof

ltmTestScore :: TrxLtmProf -> Int -> String
ltmTestScore (TrxLtmProf ltmp@(TrxProfile _ ltm) chxp@(TrxProfile _ chx)) dstart
  = intercalate "\t" . addname . maybe ["n/a"] (atstart . (+ dstart)) $ trxProfileStart ltmp
  where addname = ((BS.unpack $ trxProfileName ltmp) :)
        atstart st 
          | st < 0 = [ "n/a" ]
          | otherwise = let ltm_st = ltm U.! (st + 3)
                            ltm_ttl = U.sum ltm
                            ltm_max | (st - 1 >= 0) = (U.maxIndex . U.take 9 . U.drop (st - 1) $ ltm) - 1
                                    | otherwise     = (U.maxIndex . U.take 7                   $ ltm)
                            chx_st = chx U.! (st + 3)
                            chx_ttl = U.sum chx
                            showf :: Double -> String
                            showf x = showFFloat (Just 4) x ""
                        in [ show ltm_st, show chx_st, show ltm_ttl, show chx_ttl
                           , showf (fromIntegral ltm_st / fromIntegral ltm_ttl)
                           , showf (fromIntegral chx_st / fromIntegral chx_ttl)
                           , show ltm_max
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