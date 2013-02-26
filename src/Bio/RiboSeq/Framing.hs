module Bio.RiboSeq.Framing
       ( Terminus, TerminusIO, TerminusBase(..)
       , newTerminus, freezeTerminus
       , countAtStart, countAtEnd
         
       , Framing, FramingIO, FramingBase(..)
       , newFrame, freezeFrame
       , countInCds
       
       , lenFract, terminusPeak
       )
       where

import Control.Applicative
import Control.Monad
import Data.List
import Data.Maybe
import Data.Ord
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as VM

import qualified Bio.SamTools.Bam as Bam
import qualified Bio.SeqLoc.Location as Loc
import qualified Bio.SeqLoc.Position as Pos
import Bio.SeqLoc.Transcript

import Bio.RiboSeq.BamFile

lenFract :: Framing -> Int -> Double
lenFract fr len = (fromIntegral lenct) / (fromIntegral ttl)
  where lenct = maybe 0 V.sum $ (frprofile fr) V.!? (len - frminlen fr)
        ttl = V.sum . V.map V.sum . frprofile $ fr

terminusPeak :: (Pos.Offset, Pos.Offset) -> Terminus -> Int -> Pos.Offset
terminusPeak (minoff, maxoff) tr len = maximumBy (comparing poslen) [minoff..maxoff]
  where poslen pos = fromMaybe 0 $ do
          vlen <- (profile tr) V.!? (fromIntegral $ pos - before tr)
          vlen V.!? (len - minlen tr)

-- Count profile around a specified point in each transcript

type TerminusIO = TerminusBase (VM.IOVector Int)
type Terminus = TerminusBase (V.Vector Int)

data TerminusBase v = Tr { before, after :: !Pos.Offset
                         , minlen, maxlen :: !Int
                                   -- Outer index is (pos - before)
                                   -- Inner index is (len - minlen)
                         , profile :: !(V.Vector v)
                         }

newTerminus :: (Pos.Offset, Pos.Offset) -- ^ Nucleotides before and after the terminus
               -> (Int, Int) -- ^ Fragment length range
               -> IO TerminusIO
newTerminus (b, a) (minl, maxl) = Tr b a minl maxl <$> newProfile
  where newProfile = V.replicateM (fromIntegral $ 1 + a - b) $
                     VM.replicate (1 + maxl - minl) 0
                             
freezeTerminus :: TerminusIO -> IO Terminus                     
freezeTerminus tio = do prof <- V.mapM V.freeze . profile $ tio
                        return $ Tr { before = before tio
                                    , after = after tio
                                    , minlen = minlen tio
                                    , maxlen = maxlen tio
                                    , profile = prof
                                    }

countAtStart :: TerminusIO -> Transcript -> Bam.Bam1 -> IO ()
countAtStart mgstart trx = withCds trx $ \cdsloc ->
  onReadContig trx $ \iloc ->
  let ioff = Pos.offset . Loc.startPos $ iloc
      cdsstart = fst . Loc.bounds $ cdsloc
      len = fromIntegral . Loc.length $ iloc
  in countPosLen mgstart (ioff - cdsstart) len

countAtEnd :: TerminusIO -> Transcript -> Bam.Bam1 -> IO ()
countAtEnd mgend trx = withCds trx $ \cdsloc ->
  onReadContig trx $ \iloc ->
  let ioff = Pos.offset . Loc.startPos $ iloc
      cdsend = snd . Loc.bounds $ cdsloc
      len = fromIntegral . Loc.length $ iloc
  in countPosLen mgend (ioff - cdsend) len

countPosLen :: TerminusIO -> Pos.Offset -> Int -> IO ()
countPosLen mg pos len = when inrange $ count2D (profile mg) posidx lenidx
  where posidx = fromIntegral $ pos - before mg
        lenidx = len - minlen mg
        inrange = and [ pos >= before mg, pos <= after mg
                      , len >= minlen mg, len <= maxlen mg
                      ]

-- Frame occupancy as a function of length
        
type Framing = FramingBase (V.Vector Int)        
type FramingIO = FramingBase (VM.IOVector Int)
        
data FramingBase v = Fr { frminlen, frmaxlen :: !Int
                                    -- Outer index is (len - minlen)
                                    -- Inner index is [0..2] frame
                        , frprofile :: !(V.Vector v)
                        }

newFrame :: (Int, Int) -> IO FramingIO
newFrame (minl, maxl) = Fr minl maxl <$> newProfile
  where newProfile = V.replicateM (1 + maxl - minl) $ VM.replicate 3 0

freezeFrame :: FramingIO -> IO Framing
freezeFrame fr = do prof <- V.mapM V.freeze . frprofile $ fr
                    return $ Fr { frminlen = frminlen fr
                                , frmaxlen = frmaxlen fr
                                , frprofile = prof
                                }

countInCds :: (Pos.Offset, Pos.Offset) -> FramingIO -> Transcript -> Bam.Bam1 -> IO ()
countInCds (body5, body3) fr trx = withCds trx $ \cdsloc ->
  onReadContig trx $ \iloc ->
  let ioff = Pos.offset . Loc.startPos $ iloc
      (cdsstart, cdsend) = Loc.bounds $ cdsloc
      inCdsBody = and [ ioff >= cdsstart + body5
                      , ioff <= cdsend - body3
                      ]
      len = fromIntegral . Loc.length $ iloc          
  in when inCdsBody $ countFrame fr (ioff - cdsstart) len

countFrame :: FramingIO -> Pos.Offset -> Int -> IO ()
countFrame fr pos len = when inrange $ count2D (frprofile fr) lenidx frameidx
  where lenidx = len - frminlen fr
        frameidx = fromIntegral $ pos `mod` 3
        inrange = and [ len >= frminlen fr, len <= frmaxlen fr ]
        
count2D :: V.Vector (VM.IOVector Int) -> Int -> Int -> IO ()
count2D v i1 i2 = do x <- VM.read (v V.! i1) i2
                     VM.write (v V.! i1) i2 $! succ x
                     
