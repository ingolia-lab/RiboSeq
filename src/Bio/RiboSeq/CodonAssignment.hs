module Bio.RiboSeq.CodonAssignment ( ASites, ASiteDelta, readASiteDelta, aSiteRange, aSiteDelta, aSitePos
                                   , CodonFrame, cfCodon, cfFrame
                                   , ntOffsetToCodonFrame, codonFrameToNtOffset
                                   , ntCodonFrameFromStart, ntCodonFrameFromEnd
                                   , profileFromStart, profileFromEnd
                                   , extraBounds
                                   , trxReadContig, trxReadASite, cdsReadASite
                                   )
       where

import Control.Arrow
import Control.Applicative
import Control.Monad
import qualified Data.ByteString.Char8 as BS
import Data.Maybe

import qualified Data.Attoparsec as AP
import qualified Data.Attoparsec.Char8 as AP
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U

import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.Position as Pos
import qualified Bio.SeqLoc.SpliceLocation as SpLoc
import Bio.SeqLoc.Strand
import Bio.SeqLoc.Transcript

type ASites = [(Pos.Offset, Pos.Offset)]
type ASiteDelta = Pos.Offset -> Maybe Pos.Offset

-- | Read a mapping from fragment lengths to A site relative positions
readASiteDelta :: FilePath -> IO ASiteDelta
readASiteDelta = liftM aSiteDelta . readASite

-- | Read an association list of fragment lengths and relative A site
-- positions. 
readASite :: FilePath -> IO ASites
readASite = BS.readFile >=> either handleError return . AP.parseOnly aSiteParser
  where handleError = ioError . userError . ("Unable to parse A site file: " ++)
        
-- | Parser for a table of fragment lengths and A site relative
-- positions. Entries are lines consisting of the fragment length,
-- whitespace (e.g., a tab), and the A site offset.
aSiteParser :: AP.Parser ASites
aSiteParser = AP.many1 aSiteLine
  where aSiteLine = (,) <$> AP.decimal <*> 
                    (AP.skipSpace *> AP.signed AP.decimal <* (AP.endOfLine <|> AP.endOfInput))
                    
aSiteDelta :: ASites -> Pos.Offset -> Maybe Pos.Offset
aSiteDelta offsets = lookupLength
  where offsetVector = nothings V.// map ( fromIntegral *** Just ) offsets
          where nothings = V.replicate maxoff Nothing
                maxoff = fromIntegral . maximum . (0 :) . map (succ . fst) $ offsets
        lookupLength p = case fromIntegral p of
                              p' | p' >= V.length offsetVector -> Nothing
                                 | otherwise -> offsetVector V.! p'
                                                  
aSiteRange :: ASites -> (Pos.Offset, Pos.Offset)
aSiteRange offsets = (minimum lens, maximum lens)
  where lens = map fst offsets

-- | Determine the A site nucleotide position, which may not fall
-- precisely on a codon boundary, from the footprint alignment
-- location.
aSitePos :: (Loc.Location l) => ASiteDelta -> l -> Maybe Pos.Pos
aSitePos delta loc = delta (Loc.length loc) >>= \inpos -> Loc.posOutof (Pos.Pos inpos Plus) loc
  
data CodonFrame = CodonFrame { cfCodon, cfFrame :: !Pos.Offset }
                  deriving (Show)

-- | Take a nucleotide offset relative to the start codon and
-- determine the codon number. The codon number may be negative and
-- the function does not know the extent of the coding region. The
-- frame will be -1, 0, or 1.
ntOffsetToCodonFrame :: Pos.Offset -> CodonFrame
ntOffsetToCodonFrame = uncurry CodonFrame . second pred . flip divMod 3 . succ

-- | Invert 'ntOffsetToCodonFrame' for \frame\ of -1, 0, or 1.
codonFrameToNtOffset :: CodonFrame -> Pos.Offset
codonFrameToNtOffset (CodonFrame codon frame) = codon * 3 + frame

-- | Determine the codon position and reading frame of an A site
-- position, relative to a CDS.
ntCodonFrameFromStart :: (Loc.Location l) => l -> Pos.Pos -> Maybe CodonFrame
ntCodonFrameFromStart loc ntpos = Loc.posInto ntpos (Loc.extend (1, -1) loc) >>= toCodonFrame
  where toCodonFrame (Pos.Pos succoff Plus) = Just $! uncurry CodonFrame . second pred $! succoff `divMod` 3
        toCodonFrame (Pos.Pos _succoff Minus) = Nothing
        
-- | Determine the codon position and reading frame of an A site
-- position, relative to the /end/ of a CDS.
ntCodonFrameFromEnd :: (Stranded l, Loc.Location l) => l -> Pos.Pos -> Maybe CodonFrame
ntCodonFrameFromEnd loc ntpos = Loc.posInto ntpos (revCompl $ Loc.extend (1, -1) loc) >>= toCodonFrame
  where toCodonFrame (Pos.Pos succoff Minus) = Just $! uncurry CodonFrame . second (negate . pred) $! succoff `divMod` 3
        toCodonFrame (Pos.Pos _succoff Plus) = Nothing

-- | Convert a nucleotide occupancy profile on a transcript to a codon
-- occupancy profile on a CDS within that transcript, optionally
-- clipping it at a specified maximum length. The 
profileFromStart :: (Num a, U.Unbox a) => Maybe Int -> Loc.ContigLoc -> U.Vector a -> U.Vector a
profileFromStart mmaxlen loc ntprof = 
  let codonlen = maybe id min mmaxlen $! fromIntegral $ case Loc.length loc `divMod` 3 of
        (len, 0) -> len
        (len, _) -> len + 1
      zero = fromIntegral . Pos.offset . Loc.startPos $ loc
      atcodon c = {-# SCC "atcodon" #-}
        ( fromMaybe 0 $! ntprof U.!? (zero + 3 * c - 1) ) +
        ( fromMaybe 0 $! ntprof U.!? (zero + 3 * c) ) +
        ( fromMaybe 0 $! ntprof U.!? (zero + 3 * c + 1) )
  in U.map atcodon . U.enumFromN 0 $! codonlen
{-# SPECIALIZE profileFromStart :: Maybe Int -> Loc.ContigLoc -> U.Vector Double -> U.Vector Double #-}
{-# SPECIALIZE profileFromStart :: Maybe Int -> Loc.ContigLoc -> U.Vector Int -> U.Vector Int #-}
                             
profileFromEnd :: (Num a, U.Unbox a) => Maybe Int -> Loc.ContigLoc -> U.Vector a -> U.Vector a
profileFromEnd mmaxlen loc ntprof = 
  let codonlen = maybe id min mmaxlen $ fromIntegral $ case Loc.length loc `divMod` 3 of
        (len, 0) -> len
        (len, _) -> len + 1
      zero = fromIntegral . Pos.offset . Loc.endPos $ loc
      atcodon c = {-# SCC "atcodon_end" #-}
        ( fromMaybe 0 $! ntprof U.!? (zero - 3 * c - 3) ) +
        ( fromMaybe 0 $! ntprof U.!? (zero - 3 * c - 2) ) +
        ( fromMaybe 0 $! ntprof U.!? (zero - 3 * c - 1) )
  in U.map atcodon . U.enumFromN 0 $! codonlen
{-# SPECIALIZE profileFromEnd :: Maybe Int -> Loc.ContigLoc -> U.Vector Double -> U.Vector Double #-}
{-# SPECIALIZE profileFromEnd :: Maybe Int -> Loc.ContigLoc -> U.Vector Int -> U.Vector Int #-}

-- | Extensions of the transcript bounds to allow alignments whose A
-- site should fall within the transcript but whose full extent will
-- not.
extraBounds :: (Pos.Offset, Pos.Offset)
extraBounds = ( 100, 100 )

-- | Determine the (contiguous) transcript location covered by a
-- (possibly spliced) read location. 
--
-- In general, if the read location is incompatible with the
-- transcript structure, i.e., genomic coordinate nucleotides are
-- present in the read location that are not part of the transcript,
-- then 'Nothing' is returned. Likewise, if the read is on the 'Minus'
-- strand of the transcript then 'Nothing' is returned.
-- 
-- However, the transcript location is extended by 'extraBounds'
-- nucleotides on the start and end, in order to permit reads that
-- extend off the ends of the transcripts. This is a different
-- situation, biologically, then reads containing a splice junction
-- that is incompatible with the transcript structure.
trxReadContig :: Transcript -> SpLoc.SpliceLoc -> Maybe Loc.ContigLoc
trxReadContig trx loc = do bamiloc <- SpLoc.locInto loc exttloc
                           case Loc.toContigs bamiloc of
                             [ cloc ] | Loc.strand cloc == Plus -> Just $ Loc.slide (negate $ fst extraBounds) cloc
                             _ -> Nothing
  where (OnSeq _name tloc) = location trx
        exttloc = Loc.extend extraBounds tloc

-- | Determine the transcript location of the A site within a read.
--
-- The transcript location of the entire read is first determined by
-- 'trxReadContig'. The position of the A site within the read is
-- then determined by the 'asite' function.
--
-- If the A site itself does not lie within the transcript, then
-- 'Nothing' is returned, though portions of the overall read can
-- extend beyond the boundaries as described for 'trxReadContig'.
trxReadASite :: ASiteDelta -> Transcript -> SpLoc.SpliceLoc -> Maybe Pos.Offset
trxReadASite asite trx loc = do cloc <- trxReadContig trx loc
                                pos <- aSitePos asite cloc
                                case pos of
                                  (Pos.Pos o Plus) 
                                    | o >= 0 && o < Loc.length (unOnSeq . location $ trx) -> Just o
                                  _ -> Nothing

-- | Determine the CDS location of the A site within a read.
-- 
-- The location of the read A site on the transcript is determined as
-- described in 'trxReadASite'. The transcript coordinate is then
-- converted to a CDS coordinate. 'Nothing' is returned if the A site
-- is outside of the CDS, or if the transcript has no annotated CDS.
cdsReadASite :: ASiteDelta -> Transcript -> SpLoc.SpliceLoc -> Maybe Pos.Offset
cdsReadASite asite trx loc = cds trx >>= \cdscloc ->
  trxReadASite asite trx loc >>= \apos ->
  liftM Pos.offset $! Loc.posInto (Pos.Pos apos Plus) cdscloc
