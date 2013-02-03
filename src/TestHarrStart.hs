{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Main
  where

import Control.Applicative
import Numeric
import System.Console.GetOpt
import System.Exit
import System.FilePath
import System.IO

import qualified Data.Vector.Unboxed as U
import System.Console.CmdTheLine

import qualified Bio.SeqLoc.Bed as Bed

import Bio.RiboSeq.SVMLight
import Bio.RiboSeq.StartSVM

main :: IO ()
main = run ( testHarrStart, info)
  where info = defTI { termName = "test-harr-start"
                     , version = "0.0"
                     , termDoc = "Testing start site prediction"
                     }
        testHarrStart = test <$> binary <*> modelfile <*> bedfile <*> scorefile
        test bin modfile bed mscore = do
          modelstr <- readFile modfile
          model <- case reads modelstr of
            [(m, "")] -> return m
            _ -> do hPutStrLn stderr $! "Malformed model in " ++ show modfile
                    exitFailure
          trxs <- Bed.readBedTranscripts bed
          hPutStrLn stderr $ "Read " ++ show (length trxs) ++ " testing transcripts"
          tsc <- testHarr bin model defaultTrainPosns trxs
          putStrLn $ "TrainPosns\t" ++ displayTestCount (countPartition tsc)
          maybe (return ()) (writeScore tsc) mscore

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

modelfile :: Term String
modelfile = required $ opt Nothing $ (optInfo [ "m", "model" ])
  { optName = "MODEL", optDoc = "Model filename" }
         
bedfile :: Term String
bedfile = required $ opt Nothing $ (optInfo [ "t", "transcripts" ])
  { optName = "BED", optDoc = "Bed-format annotation filename for testing genes" }

binary :: Term String
binary = required $ opt Nothing $ (optInfo [ "b", "binary" ])
  { optName = "/PATH/TO/SVM", optDoc = "Path of directory containing SVMlight binaries" }

scorefile :: Term (Maybe String)
scorefile = value $ opt Nothing $ (optInfo [ "o", "output" ])
  { optName = "SCORE-BASE", optDoc = "Base filename for scores" }
