{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE OverloadedStrings #-}
module Bio.RiboSeq.Starts
       where

import qualified Data.ByteString.Char8 as BS
import Data.List

import qualified Data.Vector.Unboxed as U

import Bio.SeqLoc.OnSeq
import Bio.SeqLoc.Transcript

import Bio.RiboSeq.BamFile

data Start = Start { startPeak :: !Peak
                   , startTrx :: !Transcript
                   , startSequ :: !BS.ByteString
                   }
             
peakStart :: Start -> BS.ByteString
peakStart (Start pk trx sequ) = BS.intercalate "\t" $ case peakStartCodon trx sequ pk of
  [] -> noCodonFields
  [sc] -> oneCodonFields sc
  _ -> manyCodonFields
  where sharedFields = [ unSeqLabel . geneId $ trx, peakrange, peaktype, bothrange ]
        noCodonFields = sharedFields ++ [ "???" ]
        manyCodonFields = sharedFields ++ [ "***" ]
        oneCodonFields sc = sharedFields 
                            ++ [ BS.pack . show $ sc, isInBoth sc ] 
                            ++ orfAnalysis trx sequ sc
        peakrange = BS.pack $ show (ilb . peakAny $ pk) ++ "-" ++ show (iub . peakAny $ pk)
        peaktype = BS.pack $ peakTypeCode pk
        bothrange = maybe "N/A" (\i -> BS.pack $ show (ilb i) ++ "-" ++ show (iub i)) $ peakBoth pk
        isInBoth sc = BS.pack . show . maybe False (\(Interval bl bu) -> sc >= bl && sc <= bu) $ peakBoth pk

orfAnalysis trx sequ sc = [ BS.take 3 . BS.drop sc $ sequ
                          , maybe "N/A" (\cds -> BS.pack . show $ sc - cds) mcds
                          ]
  where mcds = trxVectorStart (TrxVector trx (U.empty :: U.Vector Bool))                          

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
                                                                                           
