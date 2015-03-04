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
import Data.IORef
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
                                     , fsioFailure :: !(UM.IOVector Int)
                                     , fsioTotal :: !(IORef Int)
                                     }

data FramingStats = FramingStats { fsStart, fsEnd, fsBody :: LengthFrame
                                 , fsFailure :: !(U.Vector Int)
                                 , fsTotal :: !Int
                                 } deriving (Show)

fsioNew :: (Pos.Offset, Pos.Offset) ->
           (Pos.Offset, Pos.Offset) ->
           (Int, Int) ->
           IO FramingStatsIO
fsioNew startBnds endBnds lenBnds =
  FramingStatsIO <$>
  lfioNew startBnds lenBnds <*>
  lfioNew endBnds lenBnds <*>
  lfioNew frameBnds lenBnds <*>
  UM.replicate (fromEnum (maxBound :: BamFailure)) 0 <*>
  newIORef 0
  where frameBnds = (0, 2)

fsioFreeze :: FramingStatsIO -> IO FramingStats
fsioFreeze fsio = FramingStats <$>
                  (lfioFreeze . fsioStart $ fsio) <*>
                  (lfioFreeze . fsioEnd $ fsio) <*>
                  (lfioFreeze . fsioBody $ fsio) <*>
                  (U.freeze . fsioFailure $ fsio) <*>
                  (readIORef . fsioTotal $ fsio)

fsioIncr :: FramingStatsIO -> BamFraming -> Int -> IO ()
fsioIncr fsio (Left failure) _len = let failidx = fromEnum failure
                                    in do modifyIORef' (fsioTotal fsio) succ
                                          UM.read (fsioFailure fsio) failidx >>= \n0 ->
                                            UM.write (fsioFailure fsio) failidx $! succ n0
fsioIncr fsio (Right (mstart, mend, mframe)) len = do
  modifyIORef' (fsioTotal fsio) succ
  maybe (return False) (\start -> lfioIncr (fsioStart fsio) start len) mstart
  maybe (return False) (\end -> lfioIncr (fsioEnd fsio) end len) mend
  maybe (return False) (\frame -> lfioIncr (fsioBody fsio) frame len) mframe
  return ()

data BamFailure = BamNoHit
                | BamMultiHit
                | BamFpFailure !FpFailure
                deriving (Show, Ord, Eq)

instance Bounded BamFailure where
  minBound = BamNoHit
  maxBound = BamFpFailure maxBound

instance Enum BamFailure where
  toEnum 0 = BamNoHit
  toEnum 1 = BamMultiHit
  toEnum n | n >= 2 = BamFpFailure (toEnum $ n - 2)
  fromEnum BamNoHit = 0
  fromEnum BamMultiHit = 1
  fromEnum (BamFpFailure fpf) = 2 + fromEnum fpf

type BamFraming = Either BamFailure (Maybe Pos.Offset, Maybe Pos.Offset, Maybe Pos.Offset)

-- | Find the framing of a 'Bam.Bam1' alignment of a footprint. When
-- the footprint alignment is unique -- i.e., it has a genomic
-- location as per 'Bam.refSeqLoc' and does not indicate multiple hits
-- as per 'Bam.nHits' -- then the unique genomic location is used to
-- determine the framing as per 'fpFraming'.
bamFraming :: (Pos.Offset, Pos.Offset) -> LM.SeqLocMap Transcript -> Bam.Bam1 -> BamFraming
bamFraming bodyBnds trxmap bam
  | multiHit bam = Left BamMultiHit
  | otherwise = maybe (Left BamNoHit) locFraming . Bam.refSeqLoc $ bam
  where multiHit = maybe False (> 1) . Bam.nHits
        locFraming = either (Left . BamFpFailure) Right . fpFraming bodyBnds trxmap

data FpFailure = FpNoGene
               | FpNoncodingOnly
               | FpNoncodingOverlap
               | FpMultiCoding
               | FpNoCompatible
               | FpAmbigFrame
               deriving (Show, Ord, Eq, Bounded, Enum)

type FpFraming = Either FpFailure (Maybe Pos.Offset, Maybe Pos.Offset, Maybe Pos.Offset)

groupByGene :: [Transcript] -> [[Transcript]]
groupByGene = groupBy sameGene
  where sameGene t1 t2 = geneId t1 == geneId t2

isNoncoding :: Transcript -> Bool
isNoncoding = isNothing . cds

isCoding :: Transcript -> Bool
isCoding = isJust . cds

-- | Find the framing of an @fpLoc@ specified in /genomic/ coordinates
-- relative to the 'Transcript' genome annotations that it intersects.
--
-- The full set of 'Transcript' objects whose extent overlaps with the
-- query @fpLoc@, as per 'LM.queryLocatable', is considered.
--
-- When the footprint doesn't overlap the extent of any transcript,
-- then 'FpNoGene' is returned. Otherwise, transcripts are grouped by
-- their 'geneId' values in order to analyze the gene(s) overlapping
-- the query footprint location. Multiple distinct genes yield
-- 'FpNoncodingOnly' (if all are non-coding, i.e., lack any
-- 'Transcript' with a 'cds'), 'FpNoncodingOverlap' (for a mix of
-- non-coding and coding genes), or 'FpMultiCoding' (for multiple
-- genes all coding). A single, non-coding gene hit also yields
-- 'FpNoncodingOnly'.
--
-- When all 'Transcript' objects overlapping the query @fpLoc@ share
-- the same 'geneId' and at least one of them is coding (i.e., has a
-- 'cds') then the framing is computed as per 'geneFraming'.
fpFraming :: (Pos.Offset, Pos.Offset) -> LM.SeqLocMap Transcript -> SpliceSeqLoc -> FpFraming
fpFraming bodyBnds trxmap fploc = case groupByGene $ LM.queryLocatable (Just Plus) fploc trxmap of
  [] -> Left FpNoGene
  [trxs] -> case filter isCoding trxs of
    [] -> Left FpNoncodingOnly
    codings -> geneFraming bodyBnds codings fploc
  genes -> case map (all isNoncoding) genes of
    isnc | and isnc  -> Left FpNoncodingOnly
         | or isnc   -> Left FpNoncodingOverlap
         | otherwise -> Left FpMultiCoding

-- | Find the framing of an @fpLoc@ specified in /genomic/ coordinates
-- relative to a collection of 'Transcript' objects derived from a
-- single gene.
--
-- For each individual 'Transcript' the start- and stop-relative
-- offsets are computed according to 'cdsRelIntoTranscript'. When
-- these offsets are not defined, typically because the @fpLoc@ is
-- incompatible with the splice structure of a specific 'Transcript',
-- then they are dropped from further consideration.
--
-- If no 'Transcript' has defined offsets, then 'FpNoCompatible' is
-- returned.
--
-- Reading frame positions are computed when the start-relative offset
-- is defined and no smaller than @bodyStartMin@, and the
-- stop-relative offset is defined and no larger than
-- @bodyEndMax. When more than one reading frame position is defined
-- -- i.e., there are two sets of start- and stop-relative offsets
-- falling within the specified range of offsets for the \"body\" of
-- the gene that would indicate different reading frames for @fpLoc@
-- -- then 'FpAmbigFrame' is returned.
--
-- When at least one 'Transcript' has a compatible structure with
-- @fpLoc@ and there are no conflicting reading frames within the
-- "body" of the CDS, then a successful 'FpFraming' is returned. The
-- start-relative offset is defined (i.e., not 'Nothing') when it is
-- the same across all valid transcripts, and likewise for the
-- stop-relative offset. The reading frame is defined when @fpLoc@
-- lies within the body of at least one transcript, and 'Nothing' if
-- it is not the body of any transcript; note that ambiguity in the
-- reading frame leads to an 'FpAmbigFrame' failure.
geneFraming :: (Pos.Offset, Pos.Offset) -> [Transcript] -> SpliceSeqLoc -> FpFraming
geneFraming (bodyStartMin, bodyEndMax) codings fploc =
  case mapMaybe (cdsRelIntoTranscript fploc) codings of
    [] -> Left FpNoCompatible
    termini -> let vsStart = maybeAllSame . map fst $ termini
                   vsEnd = maybeAllSame . map snd $ termini
                   frames = mapMaybe cdsRelToFrame termini
               in case group frames of
                 [] -> Right ( vsStart, vsEnd, Nothing )
                 [(fr:_)] -> Right ( vsStart, vsEnd, Just fr )
                 _ -> Left FpAmbigFrame
  where cdsRelToFrame (vsStart, vsEnd)
          | (vsStart >= bodyStartMin) && (vsEnd <= bodyEndMax) = Just $! vsStart `mod` 3
          | otherwise = Nothing
        maybeAllSame :: (Eq a) => [a] -> Maybe a
        maybeAllSame [] = Nothing
        maybeAllSame (x0:rest) = if all (== x0) rest then Just x0 else Nothing
                
-- | Find the relative position of an @fpLoc@ specified in /genomic/
-- coordinates relative to the start and stop codons of the
-- transcript, as per 'cdsRel', if these are defined. If the
-- transcript-relative coordiantes of @fpLoc@ are undefined as per
-- 'spliceSeqLocIntoContig', or the transcript lacks a 'cds', then
-- 'Nothing' is returned.
cdsRelIntoTranscript :: SpliceSeqLoc -> Transcript -> Maybe (Pos.Offset, Pos.Offset)  
cdsRelIntoTranscript fploc trx = do fpinto <- spliceSeqLocIntoContig fploc $ location trx
                                    cdsloc <- cds trx
                                    return $! cdsRel cdsloc fpinto

-- | Find a contiguous mapping of one 'SpliceSeqLoc' as a subset of an
-- outer 'SpliceSeqLoc' if it exists. This function returns 'Just' a
-- 'Loc.ContigLoc' that would satisfy @'Loc.clocOutof' into ('unOnSeq'
-- outer) == 'unOnSeq' loc@. When the reference sequence names don't
-- match, or when the query location is not a contiguous sublocation
-- of the outer location, then 'Nothing' is returned.
spliceSeqLocIntoContig (OnSeq spref sploc) (OnSeq outref outloc)
  | spref /= outref = Nothing
  | otherwise = sploc `SpLoc.locInto` outloc >>= toSingleContig
  where toSingleContig l = case Loc.toContigs l of
          [c] -> Just c
          _ -> Nothing

-- | Find the relative position of an @fpLoc@ specified in
-- /transcript-relative/ coordinates, relative to the start and stop
-- codons of a @cdsLoc@ (also defined in transcript-relative
-- coordinates). The first @Pos.Offset@ is the position relative to
-- the start codon, with '0' indicating the first nucleotide of the
-- start codon is the 5\' end of the alignment. The second
-- @Pos.Offset@ is the position relative to the stop codon, with '0'
-- indicating that the first nucleotide of the stop codon is the 5\'
-- end of the alignment.
cdsRel :: Loc.ContigLoc -> Loc.ContigLoc -> (Pos.Offset, Pos.Offset)
cdsRel cdsLoc fpLoc =
  ( Loc.offset5 fpLoc - Loc.offset5 cdsLoc
  , Loc.offset5 fpLoc - (Loc.offset5 cdsLoc + Loc.length cdsLoc - 3)
  )

