{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE OverloadedStrings #-}
module Bio.RiboSeq.Starts
       where

import Control.Monad
import qualified Data.ByteString.Char8 as BS
import Data.List
import Data.Maybe

import qualified Data.Vector.Unboxed as U

import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.Position as Pos
import qualified Bio.SeqLoc.SpliceLocation as SpLoc
import Bio.SeqLoc.OnSeq
import Bio.SeqLoc.Strand
import Bio.SeqLoc.Transcript

import Bio.RiboSeq.BamFile
import Bio.RiboSeq.Translation

data Start = Start { startPeak :: !Peak
                   , startTrx :: !Transcript
                   , startSequ :: !BS.ByteString
                   }
             
peakBed :: Start -> Transcript
peakBed (Start pk trx sequ) = Transcript { geneId = geneId trx
                                         , trxId = peakName
                                         , location = location trx
                                         , cds = peakcds
                                         }
  where peakName = toSeqLabel $ (unSeqLabel . geneId $ trx) `BS.append` "_" `BS.append` (BS.pack . ishow . peakAny $ pk)
        peakcds = case peakStartCodon trx sequ pk of
          [sc] -> let restnt = BS.drop sc sequ
                      orflen = maybe (BS.length restnt) (+ 2) . inFrameStopIdx $ restnt
                  in Just $! Loc.fromPosLen (Pos.Pos (fromIntegral sc) Plus) (fromIntegral orflen)
          _ -> Nothing

peakStart :: Start -> BS.ByteString
peakStart (Start pk trx sequ) = BS.intercalate "\t" $ case peakStartCodon trx sequ pk of
  [] -> noCodonFields
  [sc] -> oneCodonFields sc
  _ -> manyCodonFields
  where sharedFields = [ unSeqLabel . geneId $ trx, peakrange, peaktype, bothrange ]
        noCodonFields = sharedFields ++ [ "???", peakcodons ]
        manyCodonFields = sharedFields ++ [ "***", peakcodons ]
        oneCodonFields sc = sharedFields 
                            ++ [ BS.pack . show $ sc, isInBoth sc ] 
                            ++ orfAnalysis trx sequ sc
        peakrange = BS.pack . ishow . peakAny $ pk
        peaktype = BS.pack $ peakTypeCode pk
        bothrange = maybe "N/A" (BS.pack . ishow) . peakBoth $ pk
        isInBoth sc = BS.pack . show . maybe False (\(Interval bl bu) -> sc >= bl && sc <= bu) $ peakBoth pk
        peakcodons = BS.take ((ilength . peakAny $ pk) + 2) . BS.drop (ilb . peakAny $ pk) $ sequ

orfAnalysis trx sequ sc = [ BS.take 3 . BS.drop sc $ sequ
                          , mbsshow morflen
                          , mbsshow mstartoff
                          , mbsshow mendoff
                          , categ
                          , orfaa
                          ]
  where bsshow = BS.pack . show
        mbsshow = maybe "N/A" bsshow
        mcdsstart = liftM (fromIntegral . Loc.offset5) . cds $ trx
        mstartoff = liftM (sc -) mcdsstart
        mcdsend = liftM (fromIntegral . Pos.offset . Loc.endPos) . cds $ trx
        morflen = liftM (+ 2) . inFrameStopIdx . BS.drop sc $ sequ
        morfend = liftM (+ sc) morflen
        mendoff = liftM2 (-) morfend mcdsend
        orfnt = maybe (BS.drop sc sequ) (\len -> BS.take len . BS.drop sc $ sequ) morflen
        orfaa = trlOneLetter orfnt
        categ = case mstartoff of
          Nothing -> "sprc"
          Just 0 -> "canonical"
          Just startoff | startoff > 0 -> case mendoff of
            Just 0 -> "truncation"
            _      -> "internal-out-of-frame"
          Just startoff | startoff < 0 -> case (morflen, mendoff) of
            (Just len, _) | len < startoff -> "uorf"
            (_, Just 0)                    -> "extension"
            _ -> "uorf-overlapping"

peakStartCodon :: Transcript -> BS.ByteString -> Peak -> [Int]
peakStartCodon trx sequ pk | not (null augBoth) = augBoth
                           | not (null augAny)  = augAny
                           | not (null ncBoth)  = ncBoth
                           | not (null ncAny)   = ncAny
                           | otherwise = []
  where augBoth = [ i | i <- maybe [] irange $ peakBoth pk, isAug i ]
        augAny  = [ i | i <- irange $ peakAny pk, isAug i ]
        ncBoth  = [ i | i <- maybe [] irange $ peakBoth $ pk, isNC i ]
        ncAny   = [ i | i <- irange $ peakAny pk, isNC i ]
        isAug i = (BS.take 3 . BS.drop i $ sequ) == "ATG"
        isNC i = (BS.take 3 . BS.drop i $ sequ) `elem` [ "CTG", "GTG", "TTG", "ACG" ]

data Interval = Interval { ilb, iub :: !Int } deriving (Show)

irange :: Interval -> [Int]
irange (Interval ilb iub) = [ilb..iub]

ilength :: Interval -> Int
ilength (Interval ilb iub) = 1 + iub - ilb

ishow :: Interval -> String
ishow (Interval lb ub) = show lb ++ "-" ++ show ub

data Peak = Peak { peakAny :: !Interval
                 , peakHarr, peakLtm, peakBoth :: !(Maybe Interval)
                 } deriving (Show)
                            
peakTypeCode :: Peak -> String
peakTypeCode pk = case (peakBoth pk, peakHarr pk, peakLtm pk) of
  (Just _, _, _) -> "B"
  (Nothing, Just _, Nothing) -> "H"
  (Nothing, Nothing, Just _) -> "L"
  (_, _, _) -> "?"

trxPeakField :: Transcript -> Peak -> String
trxPeakField trx pk = intercalate "," . concat $ fieldses
  where fieldses = [ ival, typ, blen, start ]
        ival = [ "[" ++ show (ilb . peakAny $ pk) ++ "," ++ show (iub . peakAny $ pk) ++ "]" ]
        typ = case (peakBoth pk, peakHarr pk, peakLtm pk) of
          (Just _, _, _) -> ["B"]
          (Nothing, Just _, Nothing) -> ["H"]
          (Nothing, Nothing, Just _) -> ["L"]
          (_, _, _) -> ["?"]
        blen = case peakBoth pk of
          (Just (Interval blb bub)) -> [ show $ 1 + bub - blb ]
          _ -> []
        start = case trxVectorStart (TrxVector trx (U.empty :: U.Vector Bool)) of
          Just st | st >= (ilb . peakAny $ pk) && st <= (iub . peakAny $ pk) -> [ "%" ]
          _ -> []

peaks :: Maybe (TrxVector Bool) -> Maybe (TrxVector Bool) -> [Peak]
peaks harr ltm = maybe [] (map intervalToPeak . peakIntervals . trxVectorV) $ eitherStarts harr ltm
  where intervalToPeak iany = Peak { peakAny = iany
                                   , peakHarr = harr >>= subInterval iany
                                   , peakLtm  = ltm >>= subInterval iany
                                   , peakBoth = bothStarts harr ltm >>= subInterval iany
                                   }

eitherStarts :: Maybe (TrxVector Bool) -> Maybe (TrxVector Bool) -> Maybe (TrxVector Bool)
eitherStarts (Just tv1) (Just tv2) = Just $! TrxVector (trxVectorTrx tv1) (U.zipWith (||) (trxVectorV tv1) (trxVectorV tv2))
eitherStarts (Just tv)  Nothing    = Just tv
eitherStarts Nothing    (Just tv)  = Just tv
eitherStarts Nothing    Nothing    = Nothing

bothStarts :: Maybe (TrxVector Bool) -> Maybe (TrxVector Bool) -> Maybe (TrxVector Bool)
bothStarts (Just tv1) (Just tv2) = Just $! TrxVector (trxVectorTrx tv1) (U.zipWith (&&) (trxVectorV tv1) (trxVectorV tv2))
bothStarts _          _          = Nothing

subInterval :: Interval -> TrxVector Bool -> Maybe Interval
subInterval (Interval lb ub) tv
  = let subvec = U.take (1 + ub - lb) . U.drop lb . trxVectorV $ tv
    in case U.findIndices id subvec of
      truevec | U.null truevec -> Nothing
              | otherwise -> let slb = lb + U.head truevec
                                 sub = lb + U.last truevec
                             in Just $! Interval slb sub

peakIntervals :: U.Vector Bool -> [Interval]
peakIntervals = consecutives . U.findIndices id
  where consecutives = unfoldr splitConsec
          where splitConsec v0 | U.null v0 = Nothing
                               | otherwise = let start = U.head v0
                                                 end = U.foldl1' isConsec v0
                                                 !ival = Interval start end
                                                 rest = U.dropWhile (<= end) v0
                                             in Just (ival, rest)
                isConsec curr next | next == curr + 1 = next
                                   | otherwise = curr
                                                                                           
