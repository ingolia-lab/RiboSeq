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
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM

import qualified Bio.SamTools.Bam as Bam
import qualified Bio.SeqLoc.LocMap as LM
import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.Position as Pos
import qualified Bio.SeqLoc.SpliceLocation as SpLoc
import Bio.SeqLoc.Strand
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
                     
data LengthFrameIO = LengthFrameIO { lfioMinLen :: !Int
                                   , lfioMinPos :: !Pos.Offset
                                   , lfioProfile :: !(VM.IOVector (UM.IOVector Int))
                                   }

lfioFreeze :: LengthFrameIO -> IO (LengthFrame)
lfioFreeze lfio = LengthFrame <$>
                  pure (lfioMinLen lfio) <*>
                  pure (lfioMinPos lfio) <*>
                  ((V.freeze . lfioProfile $ lfio) >>= V.mapM U.freeze)

lfioNew :: (Pos.Offset, Pos.Offset) -> (Int, Int) -> IO LengthFrameIO
lfioNew (minpos, maxpos) (minlen, maxlen)
  = LengthFrameIO <$>
    pure minlen <*>
    pure minpos <*>
    VM.replicateM (fromIntegral $ 1 + maxpos - minpos) (UM.replicate (1 + maxlen - minlen) 0)

lfioIncr :: LengthFrameIO -> Pos.Offset -> Int -> IO Bool
lfioIncr lfio pos len = if posidx >= 0 && posidx < (VM.length . lfioProfile $ lfio)
                        then do lprof <- VM.read (lfioProfile lfio) posidx
                                if lenidx >= 0 && lenidx < (UM.length lprof)
                                  then do n0 <- UM.read lprof lenidx
                                          UM.write lprof lenidx $! succ n0
                                          return True
                                  else return False
                        else return False
  where posidx = fromIntegral $ pos - (lfioMinPos lfio)
        lenidx = len - (lfioMinLen lfio)

data LengthFrame = LengthFrame { lfMinLen :: !Int
                               , lfMinPos :: !Pos.Offset
                               , lfProfile :: !(V.Vector (U.Vector Int))
                               } deriving (Show)
                   
data FramingStatsIO = FramingStatsIO { fsioStart, fsioEnd, fsioBody :: !LengthFrameIO
                                     , fsioBodyStartMin :: Pos.Offset
                                     , fsioBodyEndMax :: Pos.Offset
                                     }

data FramingStats = FramingStats { fsStart, fsEnd, fsBody :: LengthFrame
                                 , fsBodyStartMin :: Pos.Offset
                                 , fsBodyEndMax :: Pos.Offset
                                 } deriving (Show)

fsioNew :: (Pos.Offset, Pos.Offset) ->
           (Pos.Offset, Pos.Offset) ->
           Pos.Offset -> Pos.Offset ->
           (Int, Int) ->
           IO FramingStatsIO
fsioNew startBnds endBnds bodyStartMin bodyEndMax lenBnds =
  FramingStatsIO <$>
  lfioNew startBnds lenBnds <*>
  lfioNew endBnds lenBnds <*>
  lfioNew frameBnds lenBnds <*>
  (pure bodyStartMin) <*>
  (pure bodyEndMax)
  where frameBnds = (0, 2)

fsioFreeze :: FramingStatsIO -> IO FramingStats
fsioFreeze fsio = FramingStats <$>
                  (lfioFreeze . fsioStart $ fsio) <*>
                  (lfioFreeze . fsioEnd $ fsio) <*>
                  (lfioFreeze . fsioBody $ fsio) <*>
                  (pure . fsioBodyStartMin $ fsio) <*>
                  (pure . fsioBodyEndMax $ fsio)

-- fsioIncr :: FramingStatsIO -> (Pos.Offset, Pos.Offset) -> Int -> IO (Bool, Bool, Bool)
-- fsioIncr fsio (vsStart, vsEnd) len = do
--   atStart <- lfioIncr (fsioStart fsio) vsStart len
--   atEnd <- lfioIncr (fsioEnd fsio) vsEnd len
--   inBody <- if (vsStart >= fsioBodyStartMin fsio) && (vsEnd <= fsioBodyEndMax fsio)
--             then lfioIncr (fsioBody fsio) (vsStart `mod` 3) len
--             else return False
--   return (atStart, atEnd, inBody)

-- | When the location of @bam@ (as per 'Bam.refSeqLoc') lies within
-- the location of @trx@ and on the forward strand, and @trx@ has an
-- annotated CDS, return @Just@ the position of the 5\' end of the hit
-- relative to the CDS.  The first @Pos.Offset@ is the position
-- relative to the start codon, with '0' indicating the first
-- nucleotide of the start codon is the 5\' end of the alignment. The
-- second @Pos.Offset@ is the position relative to the stop codon,
-- with '0' indicating that the first nucleotide of the stop codon is
-- the 5\' end of the alignment.

data FpFailure = FpNoGene
               | FpNoncodingOnly
               | FpNoncodingOverlap
               | FpMultiCoding
               | FpNoCompatible
               deriving (Show, Ord, Eq, Bounded, Enum)

type FpFraming = Either FpFailure (Maybe Pos.Offset, Maybe Pos.Offset, Maybe Pos.Offset)

groupByGene :: [Transcript] -> [[Transcript]]
groupByGene = groupBy sameGene
  where sameGene t1 t2 = geneId t1 == geneId t2

isNoncoding :: Transcript -> Bool
isNoncoding = isNothing . cds

isCoding :: Transcript -> Bool
isCoding = isJust . cds

maybeAllSame :: (Eq a) => [a] -> Maybe a
maybeAllSame [] = Nothing
maybeAllSame (x0:rest) = if all (== x0) rest then Just x0 else Nothing

fpFramings :: (Pos.Offset, Pos.Offset) -> LM.SeqLocMap Transcript -> SpliceSeqLoc -> FpFraming
fpFramings bodyBnds trxmap fploc = case groupByGene $ LM.queryLocatable (Just Plus) fploc trxmap of
  [] -> Left FpNoGene
  [trxs] -> case filter isCoding trxs of
    [] -> Left FpNoncodingOnly
    codings -> geneFramings bodyBnds codings fploc
  genes -> case map (all isNoncoding) genes of
    isnc | and isnc  -> Left FpNoncodingOnly
         | or isnc   -> Left FpNoncodingOverlap
         | otherwise -> Left FpMultiCoding

geneFramings :: (Pos.Offset, Pos.Offset) -> [Transcript] -> SpliceSeqLoc -> FpFraming
geneFramings (bodyStartMin, bodyEndMax) codings fploc =
  case mapMaybe (cdsRelIntoTranscript fploc) codings of
    [] -> Left FpNoCompatible
    termini -> let (vsStarts, vsEnds) = unzip termini
                   frames = mapMaybe cdsRelToFrame termini
               in Right ( maybeAllSame vsStarts
                        , maybeAllSame vsEnds
                        , maybeAllSame frames
                        )
  where cdsRelToFrame (vsStart, vsEnd)
          | (vsStart >= bodyStartMin) && (vsEnd <= bodyEndMax) = Just $! vsStart `mod` 3
          | otherwise = Nothing
                
cdsRelIntoTranscript :: SpliceSeqLoc -> Transcript -> Maybe (Pos.Offset, Pos.Offset)  
cdsRelIntoTranscript fploc trx = do fpinto <- spliceSeqLocIntoContig fploc $ location trx
                                    cdsloc <- cds trx
                                    return $! cdsRel cdsloc fpinto

spliceSeqLocIntoContig (OnSeq spref sploc) (OnSeq outref outloc)
  | spref /= outref = Nothing
  | otherwise = sploc `SpLoc.locInto` outloc >>= toSingleContig
  where toSingleContig l = case Loc.toContigs l of
          [c] -> Just c
          _ -> Nothing

cdsRel :: Loc.ContigLoc -> Loc.ContigLoc -> (Pos.Offset, Pos.Offset)
cdsRel cdsLoc fpLoc =
  ( Loc.offset5 fpLoc - Loc.offset5 cdsLoc
  , Loc.offset5 fpLoc - (Loc.offset5 cdsLoc + Loc.length cdsLoc - 3)
  )

