module Bio.RiboSeq.RrnaFrags
       where

import Control.Monad
import qualified Data.ByteString.Char8 as BS
import Data.Int
import Data.List
import Data.Maybe
import Data.Ord
import Numeric

import qualified Data.HashMap.Strict as M

import qualified Data.Vector as V
import qualified Data.Vector.Mutable as VM

import Bio.RiboSeq.RrnaStats

fragmentTable :: RrnaStats -> String
fragmentTable stats = unlines $ zipWith fragLine frags (cumulDepletion frags)
  where frags = fragments stats
        ttl = (fromIntegral $ total stats) :: Double
        showfrac n = showFFloat (Just 4) (fromIntegral n / ttl) ""
        fragLine frag cd = intercalate "\t" [ BS.unpack . refname $ frag
                                            , show (fst (fragpeak frag) + 1)
                                            , show (fst (startedge frag) + 1)
                                              ++ "-" ++ show (fst (endedge frag) + 1)
                                            , show (snd . fragpeak $ frag)
                                            , showfrac (snd . fragpeak $ frag)
                                            , showfrac cd
                                            ]

data Fragment = Fragment { refname :: !BS.ByteString
                         , fragpeak :: !(Int, Int64)
                         , startedge :: !(Int, Int64)
                         , endedge :: !(Int, Int64)
                         } deriving (Show)

fragments :: RrnaStats -> [Fragment]
fragments stats = sortFrags . concatMap seqCands . M.toList . unRrnaStats $ stats
  where seqCands (name, seqstats) = map peakfrag $ peaks stats seqstats
          where peakfrag pk = Fragment { refname = name                                    
                                       , fragpeak = pk
                                       , startedge = peakStartEdge . peakStarts seqstats $ fst pk
                                       , endedge = peakEndEdge . peakEnds seqstats $ fst pk
                                       }
        sortFrags = sortBy (comparing $ negate . snd . fragpeak)

cumulDepletion :: [Fragment] -> [Int64]
cumulDepletion = scanl1 (+) . map (snd . fragpeak)

peakEndEdge :: [(Int, Int64)] -> (Int, Int64)
peakEndEdge ends = last . filter aboveEdge $ rcumul
  where pkttl = sum . map snd $ ends
        threshold = max 1 . ceiling . (* peakEdgeCumul) . fromIntegral $ pkttl
        aboveEdge = (>= threshold) . snd
        rcumul = snd $ mapAccumR (\cumul (idx, cts) -> (cumul + cts, (idx, cumul + cts))) 0 ends

peakEnds :: V.Vector (V.Vector Int64) -> Int -> [(Int, Int64)]
peakEnds cts peak = catMaybes . V.toList . V.imap (\idx -> liftM ((,) idx)) $ ends
  where ends = V.modify setEnds $ V.replicate (V.length cts) Nothing
        setEnds e =
          V.forM_ (V.imap (,) cts) $ \(start, seqposstat) ->
          V.forM_ (V.imap (,) seqposstat) $ \(len, n) -> 
          case start + len - 1 of
            end | end < 0 || end >= VM.length e -> return ()
                | end < peak || start > peak -> return ()
                | otherwise -> do ct0 <- VM.read e end
                                  VM.write e end $! Just $! (fromMaybe 0 ct0) + n

peakStartEdge :: [(Int, Int64)] -> (Int, Int64)
peakStartEdge starts = head . filter aboveEdge $ lcumul
  where pkttl = sum . map snd $ starts
        threshold = max 1 . ceiling . (* peakEdgeCumul) . fromIntegral $ pkttl
        aboveEdge = (>= threshold) . snd
        lcumul = snd $ mapAccumL (\cumul (idx, cts) -> (cumul + cts, (idx, cumul + cts))) 0 starts

peakStarts :: V.Vector (V.Vector Int64) -> Int -> [(Int, Int64)]
peakStarts cts peak = catMaybes . V.toList . V.imap peakStartsAt $ cts
  where peakStartsAt idx lcts | idx > peak || lenToPeak >= V.length lcts = Nothing
                              | otherwise = Just (idx, V.sum . V.drop lenToPeak $ lcts)
                                where lenToPeak = peak - idx

peaks :: RrnaStats -> V.Vector (V.Vector Int64) -> [(Int, Int64)]
peaks stats = unfoldr pickpeak . coverSeqStat
  where pickpeak cover = case V.maxIndex cover of
          idx | cover V.! idx > threshold -> Just ((idx, cover V.! idx), suppress cover idx)
          _ -> Nothing
        threshold = max 1 $ ceiling $ peakMinFract * fromIntegral (total stats)
        suppress cover idx = cover V.// [ (i, 0) | i <- [minsup..maxsup] ]
          where minsup = max 0 (idx - peakHalfWidth)
                maxsup = min (V.length cover - 1) (idx + peakHalfWidth)

peakHalfWidth :: Int
peakHalfWidth = 25

peakMinFract :: Double
peakMinFract = 0.01

peakEdgeCumul :: Double
peakEdgeCumul = 0.90