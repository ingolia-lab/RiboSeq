{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Main
  where

import Control.Applicative
import System.Exit
import System.FilePath
import System.IO

import System.Console.CmdTheLine

import qualified Bio.SeqLoc.Bed as Bed

import Bio.RiboSeq.StartSVM

main :: IO ()
main = run ( trainHarrStart, info)
  where info = defTI { termName = "train-harr-start"
                     , version = "0.0"
                     , termDoc = "Training start site prediction"
                     }
        trainHarrStart = train <$> harrModel <*> bedfile <*> model
        train harrmod bed modfile = do
          trxs <- Bed.readBedTranscripts bed
          hPutStrLn stderr $ "Read " ++ show (length trxs) ++ " training transcripts"
          err <- trainHarr harrmod defaultTrainPosns trxs
          case err of 
            ExitFailure _ -> hPutStrLn stderr $ "*** SVM training failed"
            ExitSuccess -> do writeFile modfile $ show harrmod
                              hPutStrLn stderr $ "SVM trained"

instance ArgVal HarrSample where
  converter = let (pairParser, pairPrinter) = pair ','
              in ( either Left (Right . uncurry HarrSample) . pairParser
                 , (\(HarrSample bam asites) -> pairPrinter (bam, asites))
                 )

modelSvmFile :: String -> String
modelSvmFile model = (dropExtension model) ++ ".svmlight"

harrModel :: Term HarrModel
harrModel = HarrModel <$>
            pure defaultHarrAaFields <*>
            samples <*>
            (modelSvmFile <$> model) <*>
            mintotal <*>
            binary

model :: Term String
model = required $ opt Nothing $ (optInfo [ "m", "model" ])
  { optName = "MODEL", optDoc = "Model filename" }
         
samples :: Term [HarrSample]
samples = nonEmpty $ optAll [] $ (optInfo [ "h", "harr" ])
  { optName = "BAM,ASITE", optDoc = "Bam-format alignments of harringtonine treatment, with a sites file" }

bedfile :: Term String
bedfile = required $ opt Nothing $ (optInfo [ "t", "transcripts" ])
  { optName = "BED", optDoc = "Bed-format annotation filename for training genes" }

binary :: Term String
binary = required $ opt Nothing $ (optInfo [ "b", "binary" ])
  { optName = "/PATH/TO/SVM", optDoc = "Path of directory containing SVMlight binaries" }

mintotal :: Term Int
mintotal = value $ opt 50 $ (optInfo [ "z", "min-total" ])
  { optName = "MIN-COUNT", optDoc = "Minimum total reads on scored positions" }
           
