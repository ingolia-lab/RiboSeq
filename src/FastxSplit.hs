{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE Rank2Types #-}
module Main
  where

import Control.Applicative
import Control.Monad
import Control.Monad.IO.Class
import qualified Data.ByteString.Char8 as BS
import Data.Char
import qualified Data.HashMap.Strict as HM
import Data.List
import Data.Maybe
import Numeric
import System.Exit
import System.FilePath
import System.IO

import qualified Data.Vector.Unboxed as U
import System.Console.CmdTheLine

import qualified Control.Monad.Trans.Resource as R
import qualified Data.Conduit as C
import qualified Data.Conduit.Binary as CB
import qualified Data.Conduit.List as C

main :: IO ()
main = run ( fxs, info )
  where info = defTI { termName = "fastx-split"
                     , version = "0.0"
                     , termDoc = "Split FastQ file using index and random nucleotides"
                     }
        fxs = fastxSplit <$> argConf

fastxSplit :: Conf -> IO ()
fastxSplit conf = return ()

inputs :: Term [FilePath]
inputs = nonEmpty $ posAny [] (posInfo { posName = "INPUT", posDoc = "FastQ input" })

inputSource :: (MonadIO m, R.MonadResource m) => [FilePath] -> C.Producer m BS.ByteString
inputSource []    = fail "No input files"
inputSource ["-"] = CB.sourceHandle stdin
inputSource fs@(_:_) = foldl1' (*>) (map CB.sourceFile fs)

data Conf = Conf { cInputs :: ![FilePath],
                   cOutDir :: !FilePath,
                   cMinInsert :: !Int
                 } deriving (Read, Show)

argConf :: Term Conf
argConf = Conf <$>
          inputs <*>
          outDir <*>
          minInsert

outDir :: Term FilePath
outDir = required $ opt Nothing $ (optInfo [ "o", "output-dir" ])
         { optName = "OUTPUT-DIR", optDoc = "Output directory name" }

minInsert :: Term Int
minInsert = value $ opt 0 $ (optInfo [ "m", "min-insert" ])
            { optName = "MIN-INSERT", optDoc = "Minimum insert length" }

-- trxBed :: Term String
-- trxBed = required $ opt Nothing $ (optInfo [ "t", "transcripts" ])
--   { optName = "TRANSCRIPTS.BED", optDoc = "Transcript BED file" }

-- coordInput :: Term String
-- coordInput = required $ opt Nothing $ (optInfo [ "i", "input" ])
--   { optName = "INPUT.TXT", optDoc = "Input coordinates" }

-- outputFile :: Term (Maybe String)
-- outputFile = value $ opt Nothing $ (optInfo [ "o", "output" ])
--   { optName = "OUTPUT.TXT", optDoc = "Output coordinates" }

-- argConf :: Term Conf
-- argConf = Conf <$>
--           confInputType <*>
--           confStrandedness <*>
--           confCodingRelative <*>
--           confReportNoHit
--   where confInputType :: Term InputType
--         confInputType = value $ vFlag InputBed [(InputBedPE, (optInfo [ "bedpe" ]) { optName = "Bed PE (bedtools) format", optDoc = "Bed PE (bedtools) format" })]
--         confStrandedness :: Term Strandedness
--         confStrandedness = value $ vFlag FwdOnly
--                            [( FwdOnly, (optInfo [ "fwd"  ]) { optName = "Forward feature strand only", optDoc = "Forward feature strand only" }),
--                             ( RevOnly, (optInfo [ "rev"  ]) { optName = "Reverse feature strand only", optDoc = "Reverse feature strand only" }),
--                             ( Both,    (optInfo [ "both" ]) { optName = "Both feature strands",        optDoc = "Both feature strands" })]
--         confCodingRelative :: Term Bool
--         confCodingRelative = value $ vFlag Transcript 
--                              flag $ (optInfo [ "c", "cds" ]) { optName = "Coords relative to CDS start", optDoc = "Coords relative to CDS start" }
--         confReportNoHit :: Term Bool
--         confReportNoHit = value $ flag $ (optInfo [ "n", "no-hit" ]) { optName = "Report no-hit lines", optDoc = "Report no-hit lines" }

-- remapCoordLines :: (Monad m) => LM.SeqLocMap Transcript -> Conf -> C.Conduit BS.ByteString m BS.ByteString
-- remapCoordLines trxs conf = CB.lines C.=$= C.concatMapM remapCoordLine C.=$= C.map (flip BS.append "\n")
--   where remapCoordLine :: (Monad m) => BS.ByteString -> m [BS.ByteString]
--         remapCoordLine = case cInputType conf of
--           InputBed -> remapBedLine trxs conf
--           etc -> fail $ "Input type " ++ show etc ++ " not implemented"

-- remapBedLine :: (Monad m) => LM.SeqLocMap Transcript -> Conf -> BS.ByteString -> m [BS.ByteString]
-- remapBedLine trxs conf l = case BS.split '\t' l of
--   (chr:startBS:endBS:name:score:strandBS:rest)
--     -> do start <- either (const . fail $ "Bad start in " ++ show l) return $ ZP.parse decimal startBS
--           end <- either (const . fail $ "Bad end in " ++ show l) return $ ZP.parse decimal endBS
--           strand <- case strandBS of
--             "+" -> return Plus
--             "-" -> return Minus
--             _ -> fail $ "Bad strand " ++ show strandBS ++ " in " ++ show l
--           let gloc = OnSeq (toSeqLabel chr) (Loc.fromBoundsStrand start (end - 1) strand)
--           let tlocs = remapLoc trxs conf gloc
--               bedLine tsloc = BS.intercalate "\t" [ unSeqLabel . onSeqLabel $ tsloc
--                                                   , BS.pack . show . Pos.unOff . fst . Loc.bounds . unOnSeq $ tsloc
--                                                   , BS.pack . show . (+ 1) . Pos.unOff . snd . Loc.bounds . unOnSeq $ tsloc
--                                                   , name, score
--                                                   , if ((Loc.strand . unOnSeq) tsloc == Plus) then "+" else "-" ]
--           if null tlocs && cReportNoHit conf
--              then return $! [BS.intercalate "\t" $ [ "N/A", ".", ".", name, score, "." ]]
--              else return $! map bedLine tlocs

-- remapLoc :: LM.SeqLocMap Transcript -> Conf -> ContigSeqLoc -> [ContigSeqLoc]
-- remapLoc trxs conf l = mapMaybe locInto cands
--   where candList = LM.querySeqLoc l trxs
--         cands = HM.elems . HM.fromList . map (\t -> (trxId t, t)) $ candList
--         locInto t = let tsloc = location t
--                     in Loc.clocInto (unOnSeq l) (unOnSeq tsloc) >>= \tloc ->
--                     if (Loc.strand tloc == Plus && cStrandedness conf == RevOnly) ||
--                        (Loc.strand tloc == Minus && cStrandedness conf == FwdOnly)
--                        then Nothing
--                        else Just (OnSeq (trxId t) tloc)

-- locInto :: Transcript -> Conf -> ContigSeqLoc -> ContigSeqLoc
-- locinto t conf l = Loc.clocInto (unOnSeq l) (unOnSeq tsloc) >>=
--                    strandCheck >>= \tloc ->
--                    cdsAdjust (OnSeq (trxId t) tloc)
--   where tsloc = location t
--         strandCheck tloc = case cStrandedness conf of
--           FwdOnly -> if (Loc.strand tloc == Plus)  then Just tloc else Nothing
--           RevOnly -> if (Loc.strand tloc == Minus) then Just tloc else Nothing
--           Both    -> Just tloc
--         cdsAdjust = case cCodingRelative conf of
          
                                                               
                   
                   
                        
