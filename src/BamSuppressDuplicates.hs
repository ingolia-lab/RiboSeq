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
import qualified Data.HashMap.Strict as HM
import Data.Hashable
import Data.List
import Data.Maybe
import Data.Ord
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
                     , man = map P [ "Identifies and removes likely PCR duplicates. Duplicates are identified based on a nucleotide tag at the end of the read name stored in the BAM file, separated from the rest of the read name by a \"#\". Reads with no tag are not subject to deduplication. When multiple reads aligning to the same position share the same nucleotide tag, one is selected arbitrarily and written as the \"unique\" representative and, if specified, the rest are written to the file of duplicates. Optionally, the unique representative can be tagged with a \"ZD\" tag indicating the total number of duplicate reads (always 2 or more) at that position. Optionally, a table of duplicate suppression statistics can be written as a tab-separated file, tabulating the duplicate status of each distinct mapping site. In this statistics file, the first column is the total number of reads aligned to the site, the second is the number of unique reads, and the third is the count of distinct sites. Thus, \"1  1  234\" would indicate 234 distinct positions with a single unique read, \"2  2  17\" would indicate 17 distinct positions with two unique reads, and \"2  1  5\" would indicate 5 positions with a single duplicated read (2 reads total, 1 unique)." ]
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
  maybe (const $ return ()) writeStats (cStatOutput conf) $ stats
  unless (cQuiet conf) $ do
    hPutStrLn stderr $ unwords [ cBamInput conf ++ ":"
                               , "Processed", (show . nReads $ stats), "reads"
                               , "at", (show . nSites $ stats), "distinct sites."
                               ]
    hPutStrLn stderr $ unwords [ cBamInput conf ++ ":"
                               , "Suppressed", (show . nDupls $ stats), "duplicates"
                               , "at", (show . nDuplSites $ stats), "sites."
                               ]
  return ()

data SiteStat = SiteStat { ssTotal, ssUnique :: !Int } deriving (Show, Eq, Ord)

instance Hashable SiteStat where
  hashWithSalt salt (SiteStat ttl uniq) = hashWithSalt salt (ttl, uniq)

data Stats = Stats { statSites :: !(HM.HashMap SiteStat Int) } deriving (Show)

emptyStats = Stats { statSites = HM.empty }

countGroup :: Stats -> [[a]] -> Stats
countGroup s0 taggroups = let ss = SiteStat { ssTotal = sum . map length $ taggroups
                                            , ssUnique = length taggroups
                                            }
                          in Stats { statSites = HM.insertWith (+) ss 1 $ statSites s0 }

nDupls :: Stats -> Int
nDupls = HM.foldlWithKey' countDupls 0 . statSites
  where countDupls ct0 (SiteStat ttl uniq) n = ct0 + ((ttl - uniq) * n)

nDuplSites :: Stats -> Int
nDuplSites = HM.foldlWithKey' countDuplSites 0 . statSites
  where countDuplSites ct0 (SiteStat ttl uniq) n = ct0 + (if ttl > uniq then n else 0)

nSites :: Stats -> Int
nSites = HM.foldl' countSites 0 . statSites
  where countSites ct0 n = ct0 + n

nReads :: Stats -> Int
nReads = HM.foldlWithKey' countReads 0 . statSites
  where countReads ct0 (SiteStat ttl _uniq) n = ct0 + (ttl * n)

writeStats :: FilePath -> Stats -> IO ()
writeStats statfile stats = BS.writeFile statfile stattable
  where stattable = BS.unlines . map statLine . sortBy (comparing fst) . HM.toList . statSites $ stats
        statLine ((SiteStat ttl uniq), n) = BS.intercalate "\t" . map (BS.pack . show) $ [ ttl, uniq, n ]

processLocGroup :: (MonadResource m) => Conf -> Bam.OutHandle -> Maybe Bam.OutHandle -> Stats -> [Bam.Bam1] -> m Stats
processLocGroup conf hout mhdup stat0 bs = do
  let taggroups = groupByTag bs
  (uniqs, dupls) <- suppressDuplicates conf taggroups
  liftIO $ forM_ uniqs $ Bam.put1 hout
  maybe (return ()) (\hdup -> liftIO $ forM_ dupls $ Bam.put1 hdup) mhdup
  return $! countGroup stat0 taggroups

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
                 , cQuiet :: !Bool
                 } deriving (Show)

argConf :: Term Conf
argConf = Conf <$>
          argBamInput <*>
          argBamOutput <*>
          argDupOutput <*>
          argStatOutput <*>
          argAnnotate <*>
          argQuiet

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
argAnnotate = value $ flag $ (optInfo ["a", "annotate"]) { optDoc = "Annotate deduplicated reads" }

argQuiet :: Term Bool
argQuiet = value $ flag $ (optInfo ["q", "quiet"]) { optDoc = "Quiet operation" }
