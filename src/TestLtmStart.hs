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
import Bio.RiboSeq.StartSite

main :: IO ()
main = run ( testLtmStart, info)
  where info = defTI { termName = "test-ltm-start"
                     , version = "0.0"
                     , termDoc = "Testing start site prediction"
                     }
        testLtmStart = test <$> bedfile <*> ltmsamples <*> offset <*> outputfile
        test bed samps d output = do
          trxs <- Bed.readBedTranscripts bed
          hPutStrLn stderr $ "Read " ++ show (length trxs) ++ " testing transcripts"
          withFile output WriteMode $ \hout ->
            forM_ trxs $ \trx -> do
              ltmp <- readLtmProfile trx samps
              let sc = ltmTestScore ltmp d
              putStrLn sc
              hPutStrLn hout sc

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
