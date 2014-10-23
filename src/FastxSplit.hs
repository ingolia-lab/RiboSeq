{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE Rank2Types #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
module Main
  where

import Control.Applicative
import Control.Monad
import Control.Monad.IO.Class
import qualified Data.ByteString.Char8 as BS
import Data.Char
import Data.Either
import qualified Data.HashMap.Strict as HM
import Data.IORef
import Data.List
import Data.Maybe
import Numeric
import System.Directory
import System.Exit
import System.FilePath
import System.IO

import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM
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
fastxSplit conf = do createDirectory $ cOutDir conf
                     writeFile (cOutFile conf "config.txt") $ show conf
                     sspecs <- readSampleSheet conf
                     R.runResourceT $
                       do samples <- mapM (mkSample conf) sspecs
                          sfate <- sampleFate conf samples
                          inputSource (cInputs conf) C.$$ writeSamples conf sfate
                          writeStats samples

inputs :: Term [FilePath]
inputs = nonEmpty $ posAny [] (posInfo { posName = "INPUT", posDoc = "FastQ input" })

inputSource :: (MonadIO m, R.MonadResource m) => [FilePath] -> C.Producer m BS.ByteString
inputSource []    = fail "No input files"
inputSource ["-"] = CB.sourceHandle stdin
inputSource fs@(_:_) | "-" `elem` fs = fail "Cannot interleave stdin with files"
                     | otherwise = foldl1' (*>) (map CB.sourceFile fs)

writeSamples :: (MonadIO m, R.MonadResource m) => Conf -> SampleFate -> C.Sink BS.ByteString m ()
writeSamples conf sfate = toFastQ C.$= C.mapM_ (liftIO . handleSeq conf sfate)

writeStats :: (MonadIO m) => [Sample] -> m ()
writeStats _ = return ()

toFastQ :: (Monad m) => C.Conduit BS.ByteString m FastQ
toFastQ = CB.lines C.$= fqloop
  where fqloop = C.peek >>= \mnext -> case mnext of
          Nothing -> return ()
          Just _ -> C.take 4 >>= \ls -> case ls of
            [lname, lsequ, _l, lqual] -> let !fq' = FQ { fqname = BS.drop 1 lname, fqseq = lsequ, fqqual = lqual }
                                         in C.yield fq' >> fqloop
            _ -> fail $ "Partial FastQ: " ++ show ls

handleSeq :: Conf -> SampleFate -> FastQ -> IO ()          
handleSeq conf sfate fq = do
  modifyIORef' (sfCount sfate) succ
  case sequLinkers (cLinkerFormat conf) fq of
    TooShort -> tooShort
    res -> case HM.lookup (sequIndex res) (sfMap sfate) of
      Nothing -> unknown res
      Just s -> sample res s
  where tooShort = do modifyIORef' (sfShortCount sfate) succ
                      maybe (return ()) (\h -> hWriteFq h fq) $ sfShort sfate
        unknown res = do modifyIORef' (sfUnknownCount sfate) succ
                         maybe (return ()) (sample res) $ sfUnknown sfate
        sample res s = seqToSample conf s fq res

seqToSample :: Conf -> Sample -> FastQ -> LinkerRes -> IO ()
seqToSample _conf s fq0 res = hWriteFq (sHandle s) fq'
  where fq' = FQ { fqname = BS.unwords [ fqname fq0, sequIndex res, sequBarcode res ]
                 , fqseq  = sequEnce res
                 , fqqual = sequQual res
                 }

data FastQ = FQ { fqname, fqseq, fqqual :: !BS.ByteString } deriving (Show, Read)

hWriteFq :: Handle -> FastQ -> IO ()
hWriteFq h fq = BS.hPutStr h $ BS.unlines [ '@' `BS.cons` fqname fq, fqseq fq, "+", fqqual fq ]

readSampleSheet :: (MonadIO m) => Conf -> m [SampleSpec]
readSampleSheet conf = (liftIO . BS.readFile . cSampleSheet $ conf) >>= parseSampleSheet conf

parseSampleSheet :: (Monad m) => Conf -> BS.ByteString -> m [SampleSpec]
parseSampleSheet conf ss = zipWithM parseSampleSpec (BS.lines ss) ([1..] :: [Int])
  where idxlen = length (sampleIndex . cLinkerFormat $ conf)
        parseSampleSpec l lno = case BS.split ',' l of
                                 [name,idx] | BS.length idx == idxlen -> return $! SampleSpec name idx
                                            | otherwise -> fail $ "Index " ++ show idx ++ " not of length " ++ show idxlen ++ " in " ++ show l ++ " line #" ++ show lno
                                 _ -> fail $ "Bad sample spec " ++ show l ++ " line #" ++ show lno

mkSample :: (R.MonadResource m, MonadIO m) => Conf -> SampleSpec -> m Sample
mkSample conf ss = Sample ss <$>
                   allocateFileHandle (sampleOutFile conf ss) WriteMode <*>
                   (liftIO $ UM.replicate (sampleIndexIdxlen . cLinkerFormat $ conf) 0) <*>
                   (liftIO $ UM.replicate (seqBarcodeIdxlen . cLinkerFormat $ conf) 0)

sampleOutFile :: Conf -> SampleSpec -> FilePath
sampleOutFile conf ss = (cOutDir conf) </> ((BS.unpack . sName $ ss) ++ ".fastq")

tooShortFile :: Conf -> FilePath
tooShortFile conf = (cOutDir conf) </> "tooshort.fastq"

sampleFate :: (MonadIO m, R.MonadResource m) => Conf -> [Sample] -> m SampleFate
sampleFate conf samples = Fate <$>
                          indexSampleMap conf samples <*>
                          (Just <$> mkSample conf unknownSpec) <*>
                          liftIO (newIORef 0) <*>
                          (Just <$> allocateFileHandle (tooShortFile conf) WriteMode) <*>
                          liftIO (newIORef 0) <*>
                          liftIO (newIORef 0)
  where unknownSpec = SampleSpec { sName = "UnknownIndex",
                                   sIndex = BS.replicate (sampleIndexIdxlen . cLinkerFormat $ conf) 'N'
                                 }

indexSampleMap :: (MonadIO m) => Conf -> [Sample] -> m (HM.HashMap BS.ByteString Sample)
indexSampleMap conf = foldl' insertSample (return HM.empty)
  where insertSample iohm0 s = do hm' <- iohm0 >>= \hm0 -> insertPerfect hm0 s
                                  if cPerfectOnly conf
                                    then checkMismatches hm' s >> return hm'
                                    else insertMismatches hm' s
        insertPerfect hm0 s = case HM.lookup (ssIndex s) hm0 of
                               Nothing -> return $! HM.insert (ssIndex s) s hm0
                               Just t -> fail $ "Index clash " ++ show (sSpec s) ++ " and " ++ show (sSpec t) ++ " for index " ++ show (ssIndex s)
        checkMismatches hm' s = foldl' checkMismatch (return ()) (mismatches $ ssIndex s)
          where checkMismatch _ mmidx = case HM.lookup mmidx hm' of
                  Nothing -> return ()
                  Just t -> warn $ "Near-match between " ++ show (sSpec s) ++ " and " ++ show (sSpec t) ++ " for index " ++ show mmidx
        insertMismatches hm' s = foldl' insertMismatch (return hm') (mismatches $ ssIndex s)
          where insertMismatch iohmin mmidx
                  = iohmin >>=
                    \hmin -> case HM.lookup mmidx hmin of
                              Nothing -> return $! HM.insert mmidx s hmin
                              Just t -> fail $ "Index clash " ++ show (sSpec s) ++ " and " ++ show (sSpec t) ++ " for index " ++ show mmidx

data SampleFate = Fate { sfMap :: !(HM.HashMap BS.ByteString Sample)
                       , sfUnknown :: !(Maybe Sample)
                       , sfUnknownCount :: !(IORef Int)
                       , sfShort ::  !(Maybe Handle)
                       , sfShortCount :: !(IORef Int)
                       , sfCount :: !(IORef Int)
                       }

newtype Index = Index { unIndex :: Int } deriving (Show, Read, Eq, Ord, Num)

mismatches :: BS.ByteString -> [BS.ByteString]
mismatches s = concatMap mmat [0..(BS.length s - 1)]
  where mmat i = let before = BS.take i s
                     after = BS.drop (i + 1) s
                 in [ BS.concat [ before, BS.singleton ch, after] | ch <- "ACGT", ch /= BS.index s i ]

toIndex :: BS.ByteString -> Index
toIndex = BS.foldl' idxch 0
  where idxch i0 ch = let !i' = i0 * 4
                          !ich = case ch of
                                  'a' -> 0
                                  'A' -> 0
                                  'c' -> 1
                                  'C' -> 1
                                  'g' -> 2
                                  'G' -> 2
                                  't' -> 3
                                  'T' -> 3
                                  _ -> error $ "toIndex: Bad nucleotide " ++ show ch
                      in i' + ich

maxIndex :: Int -> Index
maxIndex len = (4 ^ len) - 1

lenIndex :: Int -> Int
lenIndex len = 4 ^ len

data SampleSpec = SampleSpec { sName :: !BS.ByteString, sIndex :: !BS.ByteString } deriving (Show, Read)

data Sample = Sample { sSpec :: !SampleSpec,
                       sHandle :: !Handle,
                       sIndexes :: !(UM.IOVector Int), 
                       sBarcodes :: !(UM.IOVector Int)
                     }
                
ssIndex :: Sample -> BS.ByteString
ssIndex = sIndex . sSpec

data Conf = Conf { cInputs :: ![FilePath],
                   cOutDir :: !FilePath,
                   cMinInsert :: !Int,
                   cLinkerFormat :: !LinkerFormat,
                   cSampleSheet :: !FilePath,
                   cPerfectOnly :: !Bool
                 } deriving (Read, Show)

data FormatNt = Prefix !Int | Suffix !Int deriving (Show, Read, Eq)

data LinkerFormat = LinkerFormat { prefixLength :: !Int,
                                   suffixLength :: !Int,
                                   sampleIndex :: ![FormatNt],
                                   seqBarcode :: ![FormatNt] } deriving (Show, Read, Eq)

seqBarcodeIdxlen :: LinkerFormat -> Int
seqBarcodeIdxlen = lenIndex . length . seqBarcode

sampleIndexIdxlen :: LinkerFormat -> Int
sampleIndexIdxlen = lenIndex . length . sampleIndex

data LinkerRes = TooShort
               | Res { sequEnce, sequIndex, sequBarcode, sequQual :: !BS.ByteString }
               deriving (Show, Read)

sequLinkers :: LinkerFormat -> FastQ -> LinkerRes
sequLinkers lfmt fq | BS.length (fqseq fq) < (prefixLength lfmt + suffixLength lfmt) = TooShort
                    | otherwise = let (pfx, rest) = BS.splitAt (prefixLength lfmt) $ fqseq fq
                                      (ence, sfx) = BS.splitAt (BS.length rest - suffixLength lfmt) rest
                                      qual = BS.take (BS.length ence) . BS.drop (prefixLength lfmt) $ fqqual fq
                                    in Res { sequEnce = ence
                                           , sequIndex = unlinkerize pfx sfx $ sampleIndex lfmt
                                           , sequBarcode = unlinkerize pfx sfx $ seqBarcode lfmt
                                           , sequQual = qual
                                           }
                                       
unlinkerize :: BS.ByteString -> BS.ByteString -> [FormatNt] -> BS.ByteString
unlinkerize pfx sfx = BS.pack . map unlinkerChar
  where unlinkerChar (Prefix i) = pfx `BS.index` i
        unlinkerChar (Suffix i) = sfx `BS.index` i

parseLinkerFormat :: String -> String -> LinkerFormat
parseLinkerFormat pfx sfx = LinkerFormat { prefixLength = length pfx
                                         , suffixLength = length sfx
                                         , sampleIndex = idxs
                                         , seqBarcode = bcds }
  where (idxs, bcds) = partitionEithers $ unfoldr fmt (zip [0..] pfx, zip [0..] sfx)
        fmt ([], []) = Nothing
        fmt ([], ((sidx, s):srest)) | s == 'I' = Just (Left  (Suffix sidx), ([], srest))
                                    | s == 'N' = Just (Right (Suffix sidx), ([], srest))
                                    | otherwise = error $ "Unknown linker format specifier " ++ show s ++ " in " ++ show (pfx, sfx)
        fmt (((pidx, p):prest), ss) | p == 'I' = Just (Left  (Prefix pidx), (prest, ss))
                                    | p == 'N' = Just (Right (Prefix pidx), (prest, ss))
                                    | otherwise = error $ "Unknown linker format specifier " ++ show p ++ " in " ++ show (pfx, sfx)

cOutFile :: Conf -> FilePath -> FilePath
cOutFile conf = ((cOutDir conf) </>)

argConf :: Term Conf
argConf = Conf <$>
          inputs <*>
          outDir <*>
          minInsert <*>
          linkerFormat <*>
          sampleSheet <*>
          pure False

outDir :: Term FilePath
outDir = required $ opt Nothing $ (optInfo [ "o", "output-dir" ])
         { optName = "OUTPUT-DIR", optDoc = "Output directory name" }

minInsert :: Term Int
minInsert = value $ opt 0 $ (optInfo [ "m", "min-insert" ])
            { optName = "MIN-INSERT", optDoc = "Minimum insert length" }

linkerPrefix :: Term String
linkerPrefix = value $ opt "" $ (optInfo [ "p", "prefix" ])
               { optName = "PREFIX", optDoc = "Prefix format strings" }

linkerSuffix :: Term String
linkerSuffix = value $ opt "" $ (optInfo [ "x", "suffix" ])
               { optName = "SUFFIX", optDoc = "Suffix format strings" }

linkerFormat :: Term LinkerFormat
linkerFormat = parseLinkerFormat <$> linkerPrefix <*> linkerSuffix

sampleSheet :: Term String
sampleSheet = required $ opt Nothing $ (optInfo [ "s", "sample-shet" ])
              { optName = "SAMPLESHEET.CSV", optDoc = "File name of CSV-format sample sheet" }
                   
warn :: (MonadIO m) => String -> m ()                        
warn = liftIO . hPutStrLn stderr

allocateFileHandle :: (R.MonadResource m) => FilePath -> IOMode -> m Handle
allocateFileHandle filename filemode = liftM snd $ R.allocate (openFile filename filemode) hClose
