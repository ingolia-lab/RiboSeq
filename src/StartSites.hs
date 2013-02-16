{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Main
  where

import Control.Applicative
import Control.Monad
import qualified Data.ByteString.Char8 as BS
import Data.Maybe
import Numeric
import System.Exit
import System.FilePath
import System.IO

import qualified Data.Vector.Unboxed as U
import System.Console.CmdTheLine

import qualified Bio.SamTools.FaIdx as FaIdx
import qualified Bio.SeqLoc.Bed as Bed
import Bio.SeqLoc.OnSeq
import Bio.SeqLoc.Transcript

import Bio.RiboSeq.BamFile
import Bio.RiboSeq.SVMLight
import Bio.RiboSeq.StartSVM
import Bio.RiboSeq.StartLTM
import Bio.RiboSeq.Starts

main :: IO ()
main = run ( startSites, info)
  where info = defTI { termName = "start-sites"
                     , version = "0.0"
                     , termDoc = "Testing start site prediction"
                     }
        startSites = sites <$> harrmodelfile <*> ltmmodel <*> bedfile <*> fidxfile <*> outputfile
        sites mharrfile mltm bed fidx out = do
          mharr <- readHarrModel mharrfile
          trxs <- Bed.readBedTranscripts bed
          FaIdx.withFastaIndex fidx $ \hfa -> 
            withFile out WriteMode $ \hout ->
            forM_ trxs $ processTrx mharr mltm hfa hout
              
processTrx :: (Maybe HarrModel) -> (Maybe (LtmModel, LtmSamples)) -> FaIdx.InHandle -> Handle -> Transcript -> IO ()
processTrx mharr mltm hfa hout trx = do 
  hmstart <- harrStarts mharr trx
  lmstart <- ltmStarts mltm trx
  sequ <- liftM (fromMaybe BS.empty) $ FaIdx.fetchLoc hfa (location trx)
  let pks = peaks hmstart lmstart
  forM_ (peaks hmstart lmstart) $ \pk ->
    let st = Start pk trx sequ
    in BS.hPutStrLn hout $ peakStart st
  
--  mapM_ (\f -> BS.hPutStrLn hout $ (unSeqLabel . geneId $ trx) `BS.append` "\t" `BS.append` BS.pack f) fs
--  BS.hPutStrLn hout $ jointLine trx (fromJust hmstart) (fromJust lmstart)
  
jointLine trx hstart lstart = (unSeqLabel . geneId $ trx) `BS.append` "\t" `BS.append` startchars
  where startchars = BS.pack . U.toList $ U.izipWith startchar (trxVectorV hstart) (trxVectorV lstart)
        startchar nt h l | h && l     && (Just nt) == trxVectorStart hstart = '@'
                         | h && l     && (Just nt) /= trxVectorStart hstart = '#'
                         | h && not l && (Just nt) == trxVectorStart hstart = 'H'
                         | h && not l && (Just nt) /= trxVectorStart hstart = 'h'
                         | not h && l && (Just nt) == trxVectorStart hstart = 'L'
                         | not h && l && (Just nt) /= trxVectorStart hstart = 'l'
                         | (Just nt) == trxVectorStart hstart = '!'
                         | otherwise  = '.'
  


--        sites hmod ls lm = undefined
        -- testLtmStart = test <$> ltmmodel <*> bedfile <*> ltmsamples <*> offset <> outputfile
        -- test ltmmod bed samps d output = do
        --   trxs <- Bed.readBedTranscripts bed
        --   hPutStrLn stderr $ "Read " ++ show (length trxs) ++ " testing transcripts"
        --   tsc <- testLtm ltmmod samps defaultTrainPosns trxs
        --   putStr $ concat [ show . ltmMinReads $ ltmmod, "\t", showFFloat (Just 1) (ltmMinDiff ltmmod) "", "\t"
        --                   , showFFloat (Just 1) (ltmMinRatio ltmmod) "", "\t" ]                            
        --   putStrLn $ "TrainPosns\t" ++ displayTestCount (countPartition tsc)
        --   writeScore tsc output
        --   when False $
        --     withFile output WriteMode $ \hout ->
        --     forM_ trxs $ \trx -> do
        --       ltmp <- readLtmProfile trx samps
        --       let sc = ltmTestScore ltmp d
        --       putStrLn sc
        --       hPutStrLn hout sc

          
readHarrModel :: Maybe FilePath -> IO (Maybe HarrModel)
readHarrModel mfile = maybe (return Nothing) (liftM Just . parseModel <=< readFile) $ mfile
  where parseModel mstr = case reads mstr of
          [(m, "")] -> return m
          _ -> do hPutStrLn stderr $! "Malformed model in " ++ show mfile
                  exitFailure

harrStarts :: Maybe HarrModel -> Transcript -> IO (Maybe (TrxVector Bool))
harrStarts mharr trx = maybe (return Nothing) (liftM Just . modelstarts) mharr
  where modelstarts harr = liftM (TrxVector trx) $ svmStartProfile harr trx

ltmStarts :: Maybe (LtmModel, LtmSamples) -> Transcript -> IO (Maybe (TrxVector Bool))
ltmStarts mltm trx = maybe (return Nothing) (liftM Just . modelstarts) mltm
  where modelstarts (params, samples) = do prof <- readLtmProfile trx samples
                                           return $! ltmStartProfile params prof

writeScore :: TestScore -> FilePath -> IO ()
writeScore tsc outname = do writeOne "tp" $ truePos tsc
                            writeOne "fp" $ falsePos tsc
                            writeOne "tn" $ trueNeg tsc
                            writeOne "fn" $ falseNeg tsc
  where writeOne categ = writeFile (base ++ "-" ++ categ ++ ext)
                         . unlines . map (\x -> showFFloat (Just 4) x "") . U.toList
        (base, ext) = splitExtension outname

instance ArgVal HarrSample where
  converter = let (pairParser, pairPrinter) = pair ','
              in ( either Left (Right . uncurry HarrSample) . pairParser
                 , (\(HarrSample bam asites) -> pairPrinter (bam, asites))
                 )

instance ArgVal (Maybe HarrSample) where
  converter = just

fidxfile :: Term String
fidxfile = required $ opt Nothing $ (optInfo [ "f", "fasta" ])
  { optName = "FASTA", optDoc = "Indexed fasta file" }

bedfile :: Term String
bedfile = required $ opt Nothing $ (optInfo [ "t", "transcripts" ])
  { optName = "BED", optDoc = "Bed-format annotation filename for genes" }

binary :: Term String
binary = required $ opt Nothing $ (optInfo [ "b", "binary" ])
  { optName = "/PATH/TO/SVM", optDoc = "Path of directory containing SVMlight binaries" }

ltmsamples :: Term (Maybe LtmSamples)
ltmsamples = liftM2 LtmSamples <$> ltmSample <*> chxSample
  where ltmSample = value $ opt Nothing $ (optInfo [ "l", "ltm" ]) { optName = "BAM,ASITE", optDoc = "Bam-format alignments of LTM treatment, with a sites file" }
        chxSample = value $ opt Nothing $ (optInfo [ "c", "chx" ]) { optName = "BAM,ASITE", optDoc = "Bam-format alignments of CHX treatment, with a sites file" }

outputfile :: Term String
outputfile = required $ opt Nothing $ (optInfo [ "o", "output" ])
  { optName = "OUTPUT-BASE", optDoc = "Base filename for outputs" }

harrmodelfile :: Term (Maybe String)
harrmodelfile = value $ opt Nothing $ (optInfo [ "m", "model" ])
  { optName = "MODEL", optDoc = "Harringtonine data model filename" }

ltmparams :: Term LtmModel          
ltmparams = LtmModel <$> cmdMinReads <*> cmdMinDiff <*> cmdMinRatio
  where cmdMinReads = value $ opt (ltmMinReads defaultLtmModel) $ 
                      (optInfo [ "min-reads" ]) { optName = "MIN-READS", optDoc = "Mininum LTM read count at start" }
        cmdMinDiff  = value $ opt (ltmMinDiff  defaultLtmModel) $
                      (optInfo [ "min-diff"  ]) { optName = "MIN-DIFF",  optDoc = "Minimum LTM - expected at start" }
        cmdMinRatio = value $ opt (ltmMinRatio defaultLtmModel) $
                      (optInfo [ "min-ratio" ]) { optName = "MIN-RATIO", optDoc = "Minimum LTM / expected at start" }

ltmmodel :: Term (Maybe (LtmModel, LtmSamples))
ltmmodel = (\msamp params -> liftM ((,) params) msamp) <$> ltmsamples <*> ltmparams