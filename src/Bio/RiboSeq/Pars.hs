{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE Rank2Types #-}
module Bio.RiboSeq.Pars
       where

import Control.Applicative ()
import Control.Monad
import qualified Data.ByteString.Char8 as BS
import Data.Functor.Identity
import Data.List
import Data.Maybe
import Data.Ord
import Numeric

import Control.Lens

import qualified Control.Monad.Trans.Resource as R
import qualified Data.Attoparsec.ByteString.Char8 as AP
import qualified Data.Conduit as C
import qualified Data.Conduit.Binary as C
import qualified Data.Conduit.List as C
import qualified Data.HashMap.Strict as M
import qualified Data.Vector.Unboxed as U

import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.Position as Pos

data ParsTrx = ParsTrx { ptrxName :: !SeqLabel
                       , ptrxLength :: !Pos.Offset
                       , ptrxCDS :: !(Maybe Loc.ContigLoc)
                       }
type ParsMap = M.HashMap BS.ByteString ParsTrx

type ScoreMap = M.HashMap BS.ByteString (U.Vector Double)

data TrxEntry = TE { teStart, teEnd :: !Pos.Offset }
data TrxMap = TrxMap { trxTrx, trxCds, trx3Utr, trx5Utr :: !(M.HashMap BS.ByteString TrxEntry) }

emptyTrxMap :: TrxMap
emptyTrxMap = TrxMap { trxTrx = M.empty, trxCds = M.empty, trx3Utr = M.empty, trx5Utr = M.empty }

trxFTypeMap :: BS.ByteString -> Lens' TrxMap (M.HashMap BS.ByteString TrxEntry)
trxFTypeMap ftype | ftype == "Transcript" = lensTrx
                  | ftype == "CDS" = lensCds
                  | ftype == "3UTR" = lens3Utr
                  | ftype == "5UTR" = lens5Utr
                  | otherwise = error $ "Unknown feature type " ++ show ftype
  where lensTrx = lens trxTrx (\tm trx' -> tm { trxTrx = trx' })
        lensCds = lens trxCds (\tm cds' -> tm { trxCds = cds' })
        lens3Utr = lens trx3Utr (\tm utr' -> tm { trx3Utr = utr' })
        lens5Utr = lens trx5Utr (\tm utr' -> tm { trx5Utr = utr' })
        
insertTrxLocalLine :: TrxMap -> BS.ByteString -> TrxMap
insertTrxLocalLine tm0 l = case BS.split '\t' l of
  [ geneid, ftype, startbs, endbs ] 
    -> let !te = TE { teStart = parseOffset startbs
                    , teEnd = parseOffset endbs
                    }
           parseOffset = either error id . AP.parseOnly (AP.signed AP.decimal)
           insertTE = M.insertWith duplicate geneid te
           duplicate = error $ "Duplicate for " ++ show geneid ++ ": " ++ show l
       in over (trxFTypeMap ftype) insertTE tm0
  _ -> error $ "Malformed line " ++ show l

--trxLocalLinesSink :: (Monad m) => C.GSink BS.ByteString m TrxMap
--trxLocalLinesSink = C.fold insertTrxLocalLine emptyTrxMap

readTrxLocalMap :: FilePath -> IO TrxMap
readTrxLocalMap f = R.runResourceT $ C.sourceFile f C.$$ (C.lines C.=$ C.fold insertTrxLocalLine emptyTrxMap)

toParsTrxMap :: TrxMap -> ParsMap
toParsTrxMap tm = runIdentity $ M.traverseWithKey trxToPars (trxTrx tm)
  where trxToPars :: BS.ByteString -> TrxEntry -> Identity ParsTrx
        trxToPars geneid trxte = Identity $ ParsTrx { ptrxName = toSeqLabel geneid
                                                    , ptrxLength = 1 + (teEnd trxte) - (teStart trxte)
                                                    , ptrxCDS = mcds
                                                    }
          where mcds = do cdste <- M.lookup geneid (trxCds tm)
                          Loc.clocInto (Loc.fromStartEnd (teStart cdste) (teEnd cdste)) (Loc.fromStartEnd (teStart trxte) (teEnd trxte))

readParsLocalMap :: FilePath -> IO ParsMap
readParsLocalMap = liftM toParsTrxMap . readTrxLocalMap

readParsScoreMap :: FilePath -> IO ScoreMap
readParsScoreMap f = R.runResourceT $ C.sourceFile f C.$$ (C.lines C.=$ C.fold insertScore M.empty)
  where insertScore sc0 l = case BS.split '\t' l of
          [ geneid, _len, scoresbs ] -> let !scores = U.fromList . map parseScore . BS.split ';' $ scoresbs
                                            parseScore = either error id . AP.parseOnly AP.double
                                            duplicate = error $ "Duplicate for " ++ show geneid ++ ": " ++ show l
                                        in M.insertWith duplicate geneid scores sc0
          _ -> error $ "Malformed score line " ++ show l

parsStatTable :: ParsMap -> ScoreMap -> BS.ByteString
parsStatTable trxmap scoremap = BS.unlines $ header : mapMaybe statLine geneIDs
  where geneIDs = M.keys scoremap
        statLine geneid = do trx <- M.lookup geneid trxmap
                             score <- M.lookup geneid scoremap
                             cds <- ptrxCDS trx
                             case fromIntegral . fst . Loc.bounds $ cds of
                               start | start > 0 -> Just $! statUtr5 trx score start
                                     | otherwise -> Nothing
        header = BS.intercalate "\t" [ "# YORF", "Length", "Total", "Avg", "First30", "Start30", "Max30Pos", "Max30" ]
        statUtr5 :: ParsTrx -> U.Vector Double -> Int -> BS.ByteString
        statUtr5 trx score start = BS.intercalate "\t" fields
          where fields = [ unSeqLabel . ptrxName $ trx
                         , BS.pack . show $ 1 + start
                         , showf total
                         , showf $ total / len
                         , first30
                         , start30
                         , maybe "N/A" (BS.pack . show . negate) max30off
                         , maybe "N/A" (showf . window30) max30off
                         ]
                showf x = BS.pack $ showFFloat (Just 2) x ""
                len = fromIntegral $ 1 + start
                total = U.sum . U.take start $ score
                first30 | U.length score < 30 = "N/A"
                        | otherwise = showf . U.sum . U.take 30 $ score
                start30 | start >= 15 && U.length score > start + 15 = showf . U.sum . U.take 30 . U.drop (start - 15) $ score
                        | otherwise = "N/A"
                max30off | start >= 18 = Just $ maximumBy (comparing window30) [18..start]
                         | otherwise = Nothing
                window30 i = U.sum . U.take 30 . U.drop (start - i) $ score

--                *Bio.RiboSeq.Pars> trx <- readParsLocalMap "../../../HinnebuschLab/Pars/sce_transcriptome_local.tab" 
--                *Bio.RiboSeq.Pars> sc <- readParsScoreMap "../../../HinnebuschLab/Pars/sce_Score.tab" 
--                *Bio.RiboSeq.Pars> BS.writeFile "../../../HinnebuschLab/Pars/utr5_stats.txt" $ parsStats trx sc
                
