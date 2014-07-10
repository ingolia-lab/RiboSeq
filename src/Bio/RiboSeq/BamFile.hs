module Bio.RiboSeq.BamFile 
       ( mapOverTranscripts, findTranscript
       , mapOverBams
       , onReadContig, onReadASite
       , withCds
       , transcriptNtProfile, transcriptNtLengthProfile
       , TrxVector(..), TrxProfile, trxVectorTrx, trxVectorV, trxVectorName, trxVectorStart, trxVectorLength
       , TrxProfSet(..), trxProfSetName, trxProfSetStart, totalProfile
       )
       where

import Control.Monad
import qualified Data.ByteString.Char8 as BS
import qualified Data.Iteratee as Iter
import Data.List
import System.IO

import qualified Data.Vector as V
import qualified Data.Vector.Mutable as VM
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM

import qualified Bio.SamTools.Bam as Bam
import qualified Bio.SamTools.BamIndex as BamIndex
import qualified Bio.SamTools.Iteratee as BamIter
import Bio.SeqLoc.Bed
import Bio.SeqLoc.LocRepr
import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.OnSeq
import Bio.SeqLoc.Position as Pos
import qualified Bio.SeqLoc.SpliceLocation as SpLoc
import Bio.SeqLoc.Transcript

import Bio.RiboSeq.CodonAssignment
import Bio.RiboSeq.Counting

-- | Map an 'IO' action over each transcript in a collection of @BED@
-- format files
mapOverTranscripts :: [FilePath] -> (Transcript -> IO ()) -> IO ()
mapOverTranscripts beds m = mapM_ (Iter.fileDriver bamIter) beds
  where bamIter = bedTranscriptEnum . Iter.mapM_ $ m

findTranscript :: [FilePath] -> BS.ByteString -> IO Transcript
findTranscript beds trxid = mapM (Iter.fileDriver bamIter) beds >>= unique . concat
  where bamIter = bedTranscriptEnum . Iter.joinI . Iter.filter isTrx $ Iter.stream2list 
        isTrx = (== (toSeqLabel trxid)) . trxId
        unique [trx] = return trx
        unique [] = ioError . userError $ "Could not find " ++ show trxid ++ " in " ++ show beds
        unique _ = ioError . userError $ "Multiple transcripts " ++ show trxid

-- | Map an 'IO' action over all 'Bam' alignments on a transcript. The
-- reads must map to a single contiguous 'Plus'-strand region on the
-- processed transcript.
-- 
-- A warning is printed if the transcript reference sequence is not
-- found in the alignment database.
-- 
--
-- The transcript bounds are extended by 'extraBounds' to ensure that
-- all relevant reads are captures.
mapOverBams :: BamIndex.IdxHandle -> (Bam.Bam1 -> IO ()) -> Transcript -> IO Int
mapOverBams bidx onbam trx
  = case Bam.lookupTarget (BamIndex.idxHeader bidx) tname of
      Nothing -> skipping >> return (-1)
      Just tid -> BamIter.enumIndexRegion bidx tid bnds bamiter >>= Iter.run
  where (OnSeq tlabel tloc) = location trx
        tname = unSeqLabel tlabel
        bnds = case Loc.bounds tloc of
          ( start, end ) -> ( fromIntegral $ start - fst extraBounds
                            , fromIntegral $ end + snd extraBounds )
        bamiter = Iter.foldM countBam 0
        countBam n b = onbam b >> (return $! succ n)
        skipping = hPutStrLn stderr . unwords $ [ "Skipping "
                                                , show . unSeqLabel . trxId $ trx
                                                , " on ", show tname
                                                ]
                        

-- | Takes an action that acts on a footprint alignment position
-- expressed as a contiguous location on a transcript, along with that
-- transcript, and converts it to an 'IO' action that acts on a
-- 'Bam.Bam1' alignment.
-- 
-- Read positions are converted to contiguous locations on a
-- transcript according to 'trxReadContig'.
-- 
-- Alignments that do not correspond to a single contiguous
-- 'Plus'-strand location that is entirely contained within the
-- transcript are ignored, i.e., 'return ()'
onReadContig :: (Monad m) => Transcript -> (Loc.ContigLoc -> m ()) -> Bam.Bam1 -> m ()
onReadContig trx mcloc = maybe (return ()) mcloc . (Bam.refSpLoc >=> trxReadContig trx)
{-# SPECIALIZE onReadContig :: Transcript -> (Loc.ContigLoc -> IO ()) -> Bam.Bam1 -> IO () #-}

-- | Takes an 'IO' action that acts on a nucleotide A site offset on a
-- transcript, along with that transcript and the 'ASiteDelta', and
-- converts it to an 'IO' action that acts on a 'Bam.Bam1' alignment.
--
-- Read positions are converted to A site nucleotide positions on a
-- transcript according to 'trxReadASite'.
-- 
-- Alignments that do not corresond to a single contiguous
-- 'Plus'-strand location on the transcript, or whose A site lies
-- outside of the transcript, are ignored, i.e., 'return ()'
onReadASite :: ASiteDelta -> Transcript -> (Bam.Bam1 -> Pos.Offset -> IO ()) -> Bam.Bam1 -> IO ()
onReadASite asite trx mcloc bam = maybe doUnalign doAlign $! Bam.refSpLoc bam
  where doUnalign = return ()
        doAlign refloc = case trxReadASite asite trx refloc of
          (ReadASite off) -> mcloc bam off
          Incompatible -> return ()
          Outside -> return ()
          BadLength -> return ()
        _verbose str 
          = let (OnSeq _name tloc) = location trx
                exttloc = Loc.extend extraBounds tloc
                mrepr = maybe "n/a" reprStr
            in do hPutStrLn stderr . intercalate "\t" $ 
                    [ reprStr tloc
                    , mrepr $ Bam.refSpLoc bam
                    , mrepr $ Bam.refSpLoc bam >>= flip SpLoc.locInto exttloc
                    , str
                    ]

-- | Run a monadic action using the CDS location on a 'Transcript', if
-- it has one.
withCds :: (Monad m) => Transcript -> (Loc.ContigLoc -> a -> m ()) -> a -> m ()
withCds trx mact = maybe (const $ return ()) mact . cds $ trx
{-# SPECIALIZE withCds :: Transcript -> (Loc.ContigLoc -> a -> IO ()) -> a -> IO () #-}

-- | Count reads whose A site is assigned to each nucleotide position
-- on a transcript, as determined by 'trxReadASite'
transcriptNtProfile :: ASiteDelta -> BamIndex.IdxHandle -> Transcript -> IO (U.Vector Int)
transcriptNtProfile asite bidx trx = do
  ntcts <- UM.replicate (fromIntegral . Loc.length . unOnSeq . location $ trx) 0
  let count = onReadASite asite trx $ \_bam off -> 
        let ioff = fromIntegral off
        in do ct <- UM.read ntcts ioff
              UM.write ntcts ioff $! succ ct
  _n <- mapOverBams bidx count trx
  U.freeze ntcts

-- | Count reads whose A site is assigned to each nucleotide position
-- on a transcript, as determined by 'trxReadASite', stratified by
-- length.
transcriptNtLengthProfile :: ASiteDelta -> BamIndex.IdxHandle -> Transcript -> IO (V.Vector (EnumCount Pos.Offset))
transcriptNtLengthProfile asites bidx trx = 
  let lenbnds@(lenmin, lenmax) = asdRange asites
      trxlen = fromIntegral . Loc.length . unOnSeq . location $ trx
  in do ntcts <- V.replicateM trxlen (newEnumCountIO lenbnds) >>= V.thaw
        let count = onReadASite asites trx $ \bam off ->
              let ioff = fromIntegral off
              in do lenct <- VM.read ntcts ioff
                    let len = maybe 0 Loc.length . Bam.refSpLoc $ bam
                    when (len >= lenmin && len <= lenmax) $ countEnum lenct len
        _n <- mapOverBams bidx count trx
        V.freeze ntcts >>= V.mapM freezeEnumCount
           
data TrxVector a = TrxVector !Transcript !(U.Vector a)

type TrxProfile = TrxVector Int

data TrxProfSet = TrxProfSet { transcript :: !Transcript, profiles :: ![U.Vector Int] }

trxVectorTrx :: TrxVector a -> Transcript
trxVectorTrx (TrxVector trx _prof) = trx

trxVectorV :: TrxVector a -> U.Vector a
trxVectorV (TrxVector _trx prof) = prof

trxVectorName :: TrxVector a -> BS.ByteString
trxVectorName (TrxVector trx _prof) = unSeqLabel . geneId $ trx

trxVectorStart :: TrxVector a -> Maybe Int
trxVectorStart (TrxVector trx _prof) = liftM (fromIntegral . Loc.offset5) . cds $ trx

trxVectorLength :: (U.Unbox a) => TrxVector a -> Int
trxVectorLength (TrxVector _trx prof) = U.length prof

trxProfSetName :: TrxProfSet -> BS.ByteString
trxProfSetName = unSeqLabel . geneId . transcript

trxProfSetStart :: TrxProfSet -> Maybe Int
trxProfSetStart = liftM (fromIntegral . Loc.offset5) . cds . transcript

totalProfile :: TrxProfSet -> TrxProfile
totalProfile (TrxProfSet trx profs) = TrxVector trx (foldl1 (U.zipWith (+)) profs)

