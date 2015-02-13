{-# LANGUAGE BangPatterns, OverloadedStrings, RankNTypes #-}
module Main
       where

import Control.Applicative
import Control.Monad
import Control.Monad.Trans.Resource
import Control.Monad.Reader
import Control.Monad.Error
import qualified Data.ByteString.Char8 as BS
import qualified Data.Conduit as C
import qualified Data.Conduit.List as C
import qualified Data.HashSet as HS
import Data.List
import Data.Maybe
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import Numeric
import System.Console.CmdTheLine
import System.Environment
import System.IO

import Bio.SamTools.Bam as Bam
import Bio.SamTools.Conduit as Bam
import Bio.SeqLoc.Bed
import qualified Bio.SeqLoc.Location as Loc
import qualified Bio.SeqLoc.LocMap as SLM
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.Position as Pos
import qualified Bio.SeqLoc.SpliceLocation as SpLoc
import Bio.SeqLoc.Strand
import Bio.SeqLoc.Transcript

main :: IO ()
main = run ( bamannot, info )
  where info = defTI { termName = "bam-suppress-duplicates"
                     , version = "0.0"
                     , termDoc = "Suppress duplicates in a BAM file based on barcoding nucleotides"
                     }
        bamannot = bamDupSupp <$> argConf

bamDupSupp :: Conf -> IO ()
bamDupSupp conf = do
  infile <- Bam.openBamInFile (cBamInput conf)
  let header = Bam.inHeader infile
  mhdup <- liftIO $ maybe (return Nothing) (\dupfile -> liftM Just $ Bam.openBamOutFile dupfile header) (cDupOutput conf)
  stats <- withBamOutFile (cBamOutput conf) header $ \hout ->
    runResourceT $ C.runConduit $
    Bam.sourceHandle infile C.$= C.groupBy sameLocation C.$$ C.foldM (processLocGroup conf hout mhdup) emptyStats
  maybe (return ()) (\hdup -> Bam.closeOutHandle hdup) mhdup
  return ()

data Stats = Stats { nSites, nDupls :: !Int } deriving (Show, Eq)

emptyStats = Stats { nSites = 0, nDupls = 0 }

countGroup :: Stats -> ([a], [b]) -> Stats
countGroup s0 (_uniqs, []) = s0
countGroup s0 (_uniqs, dupls@(_:_)) = Stats { nSites = 1 + nSites s0
                                            , nDupls = length dupls + nDupls s0
                                            }

processLocGroup :: (MonadResource m) => Conf -> Bam.OutHandle -> Maybe Bam.OutHandle -> Stats -> [Bam.Bam1] -> m Stats
processLocGroup conf hout mhdup stat0 bs = do
  grps@(uniqs, dupls) <- suppressDuplicates conf $ groupByTag bs
  liftIO $ forM_ uniqs $ Bam.put1 hout
  maybe (return ()) (\hdup -> liftIO $ forM_ dupls $ Bam.put1 hdup) mhdup
  return $! countGroup stat0 grps

suppressDuplicates :: (MonadIO m) => Conf -> [[Bam.Bam1]] -> m ([Bam.Bam1], [Bam.Bam1])
suppressDuplicates conf taggroups = do
  (uniqs, duplses) <- liftM unzip $ mapM suppress1 taggroups
  return (uniqs, concat duplses)
  where suppress1 [] = liftIO . ioError . userError $ "suppressDuplicates: empty tag grouping"
        suppress1 [x] = return (x, [])
        suppress1 xs@(xu:xds) = do xu' <- annotIfWanted xu (length xs)
                                   return (xu', xds)
        annotIfWanted = if cAnnotate conf
                        then annotDups
                        else \b _n -> return b
  
groupByTag :: [Bam.Bam1] -> [[Bam.Bam1]]
groupByTag = groupBy sameTag
  where sameTag b1 b2 = let t1 = tag b1
                            t2 = tag b2
                        in if BS.null t1 && BS.null t2
                           then False
                           else t1 == t2

tag :: Bam.Bam1 -> BS.ByteString
tag b = case BS.break (== '#') $ Bam.queryName b of
  (orig, tagseq) | BS.null orig -> BS.empty
                 | otherwise    -> tagseq

sameLocation :: Bam.Bam1 -> Bam.Bam1 -> Bool
sameLocation b1 b2 = and [ Bam.targetID b1 =?= Bam.targetID b2
                         , Bam.position b1 =?= Bam.position b2
                         , Bam.isReverse b1 == Bam.isReverse b2
                         , Bam.cigars b1 == Bam.cigars b2
                         ]
  where (=?=) (Just x1) (Just x2) = x1 == x2
        (=?=) (Just _x) Nothing   = False
        (=?=) Nothing   _mx2      = False

annotDups :: (MonadIO m) => Bam.Bam1 -> Int -> m (Bam.Bam1)
annotDups b ndup = liftIO $ Bam.addAuxi b dupTag ndup
  where dupTag = "ZD"

data Conf = Conf { cBamInput :: !FilePath
                 , cBamOutput :: !FilePath
                 , cDupOutput :: !(Maybe FilePath)
                 , cStatOutput :: !(Maybe FilePath)
                 , cAnnotate :: !Bool
                 } deriving (Show)

argConf :: Term Conf
argConf = Conf <$>
          argBamInput <*>
          argBamOutput <*>
          argDupOutput <*>
          argStatOutput <*>
          argAnnotate

argBamInput :: Term FilePath
argBamInput = required $ opt Nothing $ (optInfo ["i", "input"])
  { optName = "INPUT.BAM", optDoc = "BAM format input file" }

argBamOutput :: Term FilePath
argBamOutput = required $ opt Nothing $ (optInfo ["o", "output"])
  { optName = "OUTPUT.BAM", optDoc = "BAM format output file" }

argDupOutput :: Term (Maybe FilePath)
argDupOutput = value $ opt Nothing $ (optInfo ["d", "dups", "duplicates"])
  { optName = "DUPLICATES.BAM", optDoc = "BAM format file of duplicates" }

argStatOutput :: Term (Maybe FilePath)
argStatOutput = value $ opt Nothing $ (optInfo ["s", "stats", "statistics"])
  { optName = "STATS.TXT", optDoc = "Output file with duplicate statistics" }

argAnnotate :: Term Bool
argAnnotate = value $ flag $ (optInfo ["u", "text-input"]) { optDoc = "Text format input" }
