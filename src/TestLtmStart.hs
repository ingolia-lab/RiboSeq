{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Main
  where

import Control.Applicative
import Control.Monad
import qualified Data.ByteString.Char8 as BS
import Numeric
import System.Console.GetOpt
import System.Exit
import System.FilePath
import System.IO

import qualified Data.Vector.Unboxed as U
import System.Console.CmdTheLine

import qualified Bio.SeqLoc.Bed as Bed

import Bio.RiboSeq.BamFile
import Bio.RiboSeq.SVMLight
import Bio.RiboSeq.StartSVM
import Bio.RiboSeq.StartLTM

main :: IO ()
main = run ( testLtmStart, info)
  where info = defTI { termName = "test-ltm-start"
                     , version = "0.0"
                     , termDoc = "Testing start site prediction"
                     }
        testLtmStart = test <$> cmdmodel <*> bedfile <*> ltmsamples <*> offset <*> outputfile
        test ltmmod bed samps d output = do
          trxs <- Bed.readBedTranscripts bed
          hPutStrLn stderr $ "Read " ++ show (length trxs) ++ " testing transcripts"
          tsc <- testLtm ltmmod samps defaultTrainPosns trxs
          putStr $ concat [ show . ltmMinReads $ ltmmod, "\t", showFFloat (Just 1) (ltmMinDiff ltmmod) "", "\t"
                          , showFFloat (Just 1) (ltmMinRatio ltmmod) "", "\t" ]                            
          putStrLn $ "TrainPosns\t" ++ displayTestCount (countPartition tsc)
          writeScore tsc output
          when False $
            withFile output WriteMode $ \hout ->
            forM_ trxs $ \trx -> do
              ltmp <- readLtmProfile trx samps
              let sc = ltmTestScore ltmp d
              putStrLn sc
              hPutStrLn hout sc

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

bedfile :: Term String
bedfile = required $ opt Nothing $ (optInfo [ "t", "transcripts" ])
  { optName = "BED", optDoc = "Bed-format annotation filename for testing genes" }

ltmsamples :: Term LtmSamples
ltmsamples = LtmSamples <$> ltmSample <*> chxSample
  where ltmSample = lastOf $ optAll [] $ (optInfo [ "l", "ltm" ]) { optName = "BAM,ASITE", optDoc = "Bam-format alignments of LTM treatment, with a sites file" }
        chxSample = lastOf $ optAll [] $ (optInfo [ "c", "chx" ]) { optName = "BAM,ASITE", optDoc = "Bam-format alignments of CHX treatment, with a sites file" }

offset :: Term Int
offset = required $ defaultOpt (Just 0) Nothing $ (optInfo [ "d", "offset" ])
  { optName = "OFFSET", optDoc = "Offset from true start site" }

outputfile :: Term String
outputfile = required $ opt Nothing $ (optInfo [ "o", "output" ])
  { optName = "OUTPUT-BASE", optDoc = "Base filename for outputs" }

cmdmodel :: Term LtmModel          
cmdmodel = LtmModel <$> cmdMinReads <*> cmdMinDiff <*> cmdMinRatio
  where cmdMinReads = value $ opt (ltmMinReads defaultLtmModel) $ 
                      (optInfo [ "min-reads" ]) { optName = "MIN-READS", optDoc = "Mininum LTM read count at start" }
        cmdMinDiff  = value $ opt (ltmMinDiff  defaultLtmModel) $
                      (optInfo [ "min-diff"  ]) { optName = "MIN-DIFF",  optDoc = "Minimum LTM - expected at start" }
        cmdMinRatio = value $ opt (ltmMinRatio defaultLtmModel) $
                      (optInfo [ "min-ratio" ]) { optName = "MIN-RATIO", optDoc = "Minimum LTM / expected at start" }
