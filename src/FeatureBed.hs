{-# LANGUAGE BangPatterns #-}
module Main
  where

import Control.Monad
import qualified Data.ByteString.Char8 as BS
import qualified Data.ByteString.Lazy.Char8 as LBS
import qualified Data.HashMap.Strict as HM
import qualified Data.HashSet as HS
import Data.List
import Data.Maybe
import System.Environment
import System.IO

import qualified Data.Iteratee.IO as IterIO
import qualified Data.Iteratee.ListLike as IterLL

import qualified Bio.SeqLoc.Bed as Bed
import Bio.SeqLoc.LocRepr
import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.Position as Pos
import qualified Bio.SeqLoc.SpliceLocation as SpLoc
import Bio.SeqLoc.Strand
import Bio.SeqLoc.Transcript

type TrxMap = HM.HashMap BS.ByteString [Transcript]

main :: IO ()
main = getArgs >>= mainWithArgs
  where mainWithArgs [ trxbed, canonkgids, outbase ] = doFeatureBed trxbed canonkgids outbase
        mainWithArgs _ = do prog <- getProgName
                            hPutStrLn stderr $ unwords [ "USAGE:", prog, "<ALL-BED>", "<CANON-ID-LIST>", "<OUTBASE>" ]

doFeatureBed :: FilePath -> FilePath -> FilePath -> IO ()
doFeatureBed trxbed canonkgids outbase = do
  trxs <- readBedRefMap trxbed
  canonids <- liftM (HS.fromList . BS.lines) . BS.readFile $ canonkgids
  let iscanon trx = HS.member (unSeqLabel . geneId $ trx) canonids
      canons = filter iscanon $ allTrxs trxs
  withFile (outbase ++ "_log.txt") WriteMode $ \hlog -> do
    withFile (outbase ++ "_cds.bed") WriteMode $ \hcds ->
      forM_ canons $ handleCds hcds hlog trxs
    withFile (outbase ++ "_stop.bed") WriteMode $ \hstop ->
      forM_ canons $ handleStop hstop hlog trxs
    withFile (outbase ++ "_utr5.bed") WriteMode $ \hutr5 ->
      forM_ canons $ handleUtr5 hutr5 hlog trxs
    withFile (outbase ++ "_utr3.bed") WriteMode $ \hutr3 ->
      forM_ canons $ handleUtr3 hutr3 hlog trxs
    withFile (outbase ++ "_intron.bed") WriteMode $ \hintron ->
      forM_ canons $ handleIntrons hintron hlog trxs
            
handleCds :: Handle -> Handle -> TrxMap -> Transcript -> IO ()            
handleCds hout hlog _trxs trx = maybe noCds doCds $ cdsLocation trx
  where doCds cdsloc = BS.hPutStrLn hout $ Bed.transcriptToBedStd cdstrx
          where cdstrx = Transcript { geneId = SeqLabel . LBS.pack $ cdsname
                                    , trxId = SeqLabel . LBS.pack $ cdsname
                                    , location = cdsloc
                                    , cds = Nothing 
                                    }
        noCds = do hPutStrLn hlog $ intercalate "\t" [ cdsname, "no cds" ]
        cdsname = (BS.unpack . unSeqLabel . geneId $ trx) ++ "_cds"
          
handleStop :: Handle -> Handle -> TrxMap -> Transcript -> IO ()
handleStop hout hlog _trxs trx = maybe noCds doStop $ cds trx >>= atStop >>= flip Loc.clocOutof trxsploc
  where noCds = do hPutStrLn hlog $ intercalate "\t" [ trxname, "no CDS" ]
        stopUpstream = 27
        stopDownstream = 25
        atStop cdsloc = let (Pos.Pos endoff _) = Loc.endPos cdsloc
                        in return $! Loc.fromStartEnd (endoff - stopUpstream) (endoff + stopDownstream)
        (OnSeq trxref trxsploc) = location trx
        trxname = BS.unpack . unSeqLabel . geneId $ trx
        doStop stoploc = 
          let stopName = trxname ++ "_stop"
              stopTrx = Transcript { geneId = SeqLabel . LBS.pack $ stopName
                                   , trxId = SeqLabel . LBS.pack $ stopName
                                   , location = (OnSeq trxref stoploc)
                                   , cds = Nothing
                                   }
          in BS.hPutStrLn hout $ Bed.transcriptToBedStd stopTrx

handleUtr3 :: Handle -> Handle -> TrxMap -> Transcript -> IO ()
handleUtr3 hout hlog trxs trx = maybe noUtr doUtr $ utr3 trx >>= shrinkUtr3 >>= flip Loc.clocOutof trxsploc
  where noUtr = do hPutStrLn hlog $ intercalate "\t" [ trxname, "no 3'UTR" ]
        shrinkOffset = 3
        shrinkUtr3 cloc = case Loc.length cloc - shrinkOffset of
          len' | len' <= 0 -> Nothing
               | otherwise -> Loc.clocOutof (Loc.fromPosLen (Pos.Pos shrinkOffset Plus) len') cloc
        (OnSeq trxref trxsploc) = location trx
        trxname = BS.unpack . unSeqLabel . geneId $ trx
        doUtr utrloc = 
          let utrSeqloc = OnSeq trxref utrloc
              utrName = trxname ++ "_3utr"
              utrTrx = Transcript { geneId = SeqLabel . LBS.pack $ utrName
                                  , trxId = SeqLabel . LBS.pack $ utrName
                                  , location = (OnSeq trxref utrloc)
                                  , cds = Nothing
                                  }
          in case find (maybe False (overlaps utrSeqloc) . cdsLocation) $ lookupTrxs trxref trxs of
            Just _overcds -> hPutStrLn hlog $ intercalate "\t" $
                             [ utrName, reprStr utrSeqloc, trxname ]
            Nothing -> BS.hPutStrLn hout $ Bed.transcriptToBedStd utrTrx

handleUtr5 :: Handle -> Handle -> TrxMap -> Transcript -> IO ()
handleUtr5 hout hlog trxs trx = maybe noUtr doUtr $ utr5 trx >>= shrinkUtr5 >>= flip Loc.clocOutof trxsploc
  where noUtr = do hPutStrLn hlog $ intercalate "\t" [ trxname, "no 5'UTR" ]
        shrinkOffset = 3
        shrinkUtr5 cloc = case Loc.length cloc - shrinkOffset of
          len' | len' <= 0 -> Nothing
               | otherwise -> Loc.clocOutof (Loc.fromPosLen (Pos.Pos 0 Plus) len') cloc
        (OnSeq trxref trxsploc) = location trx
        trxname = BS.unpack . unSeqLabel . geneId $ trx
        doUtr utrloc = 
          let utrSeqloc = OnSeq trxref utrloc
              utrName = trxname ++ "_5utr"
              utrTrx = Transcript { geneId = SeqLabel . LBS.pack $ utrName
                                  , trxId = SeqLabel . LBS.pack $ utrName
                                  , location = (OnSeq trxref utrloc)
                                  , cds = Nothing
                                  }
          in case find (maybe False (overlaps utrSeqloc) . cdsLocation) $ lookupTrxs trxref trxs of
            Nothing -> BS.hPutStrLn hout $ Bed.transcriptToBedStd utrTrx
            Just overcds -> do hPutStrLn hlog utrName
                               mapM_ (hPutStrLn hlog . ("\t" ++)) $
                                 [ reprStr utrSeqloc
                                 , maybe "N/A" reprStr $ utr5 trx
                                 , reprStr $ location trx
                                 , BS.unpack . unSeqLabel . geneId $ overcds
                                 , maybe "N/A" reprStr $ cdsLocation overcds 
                                 ]

handleIntrons :: Handle -> Handle -> TrxMap -> Transcript -> IO ()
handleIntrons hout hlog trxs trx = zipWithM_ handleIntron introns ([1..] :: [Int])
  where introns = map intron . junctions $ trxsploc
        (OnSeq trxref trxsploc) = location trx
        handleIntron intrLoc intrNo =
          let intrSeqloc = OnSeq trxref (fromJust $ SpLoc.fromContigs [ intrLoc ])
              intrName = (BS.unpack . unSeqLabel . geneId $ trx) ++ "_intron" ++ show intrNo
              intrTrx = Transcript { geneId = SeqLabel . LBS.pack $ intrName
                                   , trxId = SeqLabel . LBS.pack $ intrName
                                   , location = intrSeqloc
                                   , cds = Nothing
                                   }
          in case find (overlaps intrSeqloc . location) $ lookupTrxs trxref trxs of
            Just overtrx -> hPutStrLn hlog $ intercalate "\t" $ 
                            [ intrName, reprStr intrSeqloc, BS.unpack . unSeqLabel . geneId $ overtrx ]
            Nothing -> BS.hPutStrLn hout $ Bed.transcriptToBedStd intrTrx
          
overlaps :: SpliceSeqLoc -> SpliceSeqLoc -> Bool
overlaps (OnSeq ref loc) (OnSeq qref qsploc)
  | qref /= ref = False
  | lb > qub = False
  | ub < qlb = False
  | otherwise = Loc.overlaps loc qsploc
  where (!qlb, !qub) = Loc.bounds qsploc
        (!lb, !ub) = Loc.bounds loc
        
readBedRefMap :: FilePath -> IO (HM.HashMap BS.ByteString [Transcript])
readBedRefMap = IterIO.fileDriver bedIter
  where bedIter = Bed.bedTranscriptEnum $ IterLL.foldl' insertTrx HM.empty
        insertTrx m0 trx = let ref = unSeqLabel . onSeqLabel . location $ trx
                               !l' = trx : (HM.lookupDefault [] ref m0)
                           in HM.insert ref l' m0

lookupTrxs :: SeqLabel -> TrxMap -> [Transcript]
lookupTrxs reflabel = HM.lookupDefault [] (unSeqLabel reflabel)

allTrxs :: TrxMap -> [Transcript]
allTrxs = concat . HM.elems