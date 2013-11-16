{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Main
  where

import Control.Applicative
import Control.Monad
import qualified Data.ByteString.Char8 as BS
import Data.Char
import Data.Maybe
import Numeric
import System.Exit
import System.FilePath
import System.IO

import qualified Data.Vector.Unboxed as U
import System.Console.CmdTheLine

import Bio.RiboSeq.Pars

main :: IO ()
main = run ( parsStats, info )
  where info = defTI { termName = "pars-stats"
                     , version = "0.0"
                     , termDoc = "Table of 5'UTR PARS statistics"
                     }
        parsStats = stats <$> localTrxFile <*> scoreFile <*> outputFile
        stats trxF scoreF outputF = do
          trx <- readParsLocalMap trxF
          score <- readParsScoreMap scoreF
          BS.writeFile outputF $ parsStatTable trx score

localTrxFile :: Term String
localTrxFile = required $ opt Nothing $ (optInfo [ "l", "local-trx" ])
  { optName = "TRANSCRIPT-LOCAL.TXT", optDoc = "PARS local transcript coordinate table" }
               
scoreFile :: Term String
scoreFile = required $ opt Nothing $ (optInfo [ "s", "score" ])
  { optName = "SCORE.TXT", optDoc = "PARS score table" }

outputFile :: Term String
outputFile = required $ opt Nothing $ (optInfo [ "o", "output" ])
  { optName = "OUTPUT.txt", optDoc = "Output filename" }