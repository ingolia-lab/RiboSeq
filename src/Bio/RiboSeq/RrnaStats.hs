{-# LANGUAGE RankNTypes #-}
module Bio.RiboSeq.RrnaStats
       where

import Control.Applicative
import Control.Monad
import Control.Monad.ST
import qualified Data.ByteString.Char8 as BS
import Data.IORef
import Data.Int
import Data.List
import qualified Data.HashMap.Strict as M
import Data.Maybe
import qualified Data.Traversable as T
import Numeric

import qualified Data.Vector as V
import qualified Data.Vector.Mutable as VM

import qualified Bio.SamTools.Bam as Bam

-- | Count of rRNA alignments
data RrnaStats = RrnaStats { unRrnaStats :: !(M.HashMap BS.ByteString (V.Vector (V.Vector Int64)))
                             -- ^ 'HashMap' of reference sequence names,
                             -- outer 'V.Vector' indexed by position and
                             -- inner 'V.Vector' indexed by read length
                           , rrnaUnmap :: !Int64 -- ^ Count of unaligned reads
                           }

-- | IO-mutable version of 'RrnaStats'
data RrnaStatsIO = RrnaStatsIO { unRrnaStatsIO :: !(M.HashMap BS.ByteString (V.Vector (VM.IOVector Int64))) 
                             -- ^ 'HashMap' of reference sequence names,
                             -- outer 'V.Vector' indexed by position and
                             -- inner 'VM.Vector' indexed by read length
                               , rrnaUnmapIO   :: !(IORef Int64) -- ^ Count of unaligned reads
                               }
          
-- | Copy 'RrnaStatsIO' into a pure immutable set of counts
freezeRrnaStatsIO :: RrnaStatsIO -> IO RrnaStats
freezeRrnaStatsIO (RrnaStatsIO statsio unmapio) = RrnaStats <$> 
                                                  T.mapM freezeSeqStat statsio <*>
                                                  readIORef unmapio
  where freezeSeqStat :: V.Vector (VM.IOVector Int64) -> IO (V.Vector (V.Vector Int64))
        freezeSeqStat = V.mapM V.freeze

-- | Create a new 'RrnaStatsIO' for rRNA sequences extracted from 'Bam.HeaderSeq'
mkRrnaStatsIO :: Int -> [Bam.HeaderSeq] -> IO RrnaStatsIO
mkRrnaStatsIO maxlen hseqs = RrnaStatsIO <$> foldM insertSeqStat M.empty hseqs <*> newIORef 0
  where insertSeqStat m0 hseq = do seqstat <- mkSeqStat . fromIntegral . Bam.len $ hseq
                                   return $! M.insert (Bam.name hseq) seqstat m0
        mkSeqStat len = V.replicateM len mkSeqPosStat
        mkSeqPosStat = VM.replicate (maxlen + 1) 0
        
-- | Tally one single 'Bam.Bam1' alignment.
countAlign :: RrnaStatsIO -> Bam.Bam1 -> IO ()
countAlign stats ba | Bam.isReverse ba = incrUnmap
                    | Bam.isUnmap ba = incrUnmap
                    | otherwise = maybe incrUnmap (incr seqPosStat . fromIntegral) (Bam.queryLength ba)
  where incr vm idx = do n0 <- VM.read vm idx
                         VM.write vm idx $! n0 + 1
        seqStat = fromMaybe noSeq $ Bam.targetName ba >>= flip M.lookup (unRrnaStatsIO stats)
          where noSeq = error $ "Could not find " ++ show (Bam.targetName ba)
        seqPosStat = fromMaybe noPos $ Bam.position ba >>= (V.!?) seqStat . fromIntegral
          where noPos = error $ "Out of bounds: " ++ show (Bam.targetName ba, Bam.position ba)
        incrUnmap = {-# SCC "incrUnmap" #-}
          do ct <- readIORef (rrnaUnmapIO stats)
             writeIORef (rrnaUnmapIO stats) $! succ ct

total :: RrnaStats -> Int64
total stats = totalmap stats + rrnaUnmap stats

totalmap :: RrnaStats -> Int64
totalmap = fromIntegral . sum . map (V.sum . V.map V.sum . snd) . M.toList . unRrnaStats

-- | Collapse a position / length count vector to a starting position
-- vector
startSeqStat :: V.Vector (V.Vector Int64) -> V.Vector Int64
startSeqStat = V.map V.sum

-- | Collapse a position / length count vector to an ending position
-- vector.
endSeqStat :: V.Vector (V.Vector Int64) -> V.Vector Int64
endSeqStat seqstat = V.modify modifyEnd $ V.replicate (V.length seqstat) 0
  where modifyEnd e = V.forM_ (V.imap (,) seqstat) $ \(start, seqposstat) ->
          V.forM_ (V.imap (,) seqposstat) $ \(len, n) -> 
          case start + len - 1 of
            end | end < 0 || end >= VM.length e -> return ()
                | otherwise -> do ct0 <- VM.read e end
                                  VM.write e end $! ct0 + n

coverSeqStat :: V.Vector (V.Vector Int64) -> V.Vector Int64
coverSeqStat seqstat = V.modify modifyCover $ V.replicate (V.length seqstat) 0
  where modifyCover c = V.forM_ (V.imap (,) seqstat) $ \(start, seqposstat) ->
          V.forM_ (V.imap (,) seqposstat) $ \(len, n) -> 
          forM_ [0..(len-1)] $ \dpos ->
          case start + dpos of
            pos | pos < 0 || pos >= VM.length c -> return ()
                | otherwise -> do ct0 <- VM.read c pos
                                  VM.write c pos $! ct0 + n

countTable :: RrnaStats -> String
countTable stats = unlines [ showStat "Input" (total stats)
                           , showStat "rRNA" (totalmap stats)
                           , showStat "Kept"  (rrnaUnmap stats)
                           ]
  where ttl = (fromIntegral $ total stats) :: Double
        showStat name ct = intercalate "\t" [ name, show ct, showFFloat (Just 3) (fromIntegral ct / ttl) "" ]

posTable :: RrnaStats -> String
posTable stats = unlines . concatMap seqTable $ assocs
  where assocs = M.toList . unRrnaStats $ stats
        ttl = (fromIntegral $ total stats) :: Double
        seqTable (name, seqstat) = map seqPosLine [0..maxpos]
          where maxpos = (V.length seqstat) - 1
                starts = startSeqStat seqstat
                ends = endSeqStat seqstat
                covers = coverSeqStat seqstat
                showfrac n = showFFloat (Just 6) (fromIntegral n / ttl) ""
                seqPosLine pos = unfields [ BS.unpack name, show pos
                                          , show nstart, showfrac nstart
                                          , show nend, showfrac nend
                                          , show ncover, showfrac ncover
                                          ]
                  where nstart = fromMaybe 0 $ starts V.!? pos
                        nend = fromMaybe 0 $ ends V.!? pos
                        ncover = fromMaybe 0 $ covers V.!? pos

posLenTable :: (Int, Int) -> RrnaStats -> String
posLenTable (minlen, maxlen) stats = unlines . concatMap seqTable $ assocs
  where assocs = M.toList . unRrnaStats $ stats
        seqTable (name, seqstat) = map seqPosLine [0..maxpos]
          where maxpos = V.length seqstat - 1
                seqPosLine pos = unfields $ [ BS.unpack name, show pos ] ++ seqPosLenFields
                  where seqPosLenFields = map seqPosLenField [minlen..maxlen]
                        seqPosLenField len = maybe "0" show $ seqstat V.!? pos >>= flip (V.!?) len

unfields :: [String] -> String
unfields = intercalate "\t"
        
