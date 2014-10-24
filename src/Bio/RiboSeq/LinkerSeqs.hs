{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE Rank2Types #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
module Bio.RiboSeq.LinkerSeqs
       where

import qualified Data.ByteString.Char8 as BS
import Data.Either
import Data.List
import qualified Data.Vector.Unboxed.Mutable as UM
import System.IO

--
data FastQ = FQ { fqname, fqseq, fqqual :: !BS.ByteString } deriving (Show, Read)

hWriteFq :: Handle -> FastQ -> IO ()
hWriteFq h fq = BS.hPutStr h $ BS.unlines [ '@' `BS.cons` fqname fq, fqseq fq, "+", fqqual fq ]

--
data FormatNt = Prefix !Int | Suffix !Int deriving (Show, Read, Eq)

data LinkerFormat = LinkerFormat { prefixLength :: !Int,
                                   suffixLength :: !Int,
                                   sampleIndex :: ![FormatNt],
                                   seqBarcode :: ![FormatNt] } deriving (Show, Read, Eq)

seqBarcodeIdxlen :: LinkerFormat -> Int
seqBarcodeIdxlen = lenIndex . length . seqBarcode

sampleIndexIdxlen :: LinkerFormat -> Int
sampleIndexIdxlen = lenIndex . length . sampleIndex

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

--unparseLinkerFormat :: LinkerFormat -> (String, String)
--unparseLinkerFormat lf = 

--
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

--
newtype Index = Index { unIndex :: Int } deriving (Show, Read, Eq, Ord, Num, Integral, Real, Enum)

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

fromIndex :: Int -> Index -> BS.ByteString
fromIndex len idx = BS.reverse $! BS.unfoldr go (len, idx)
  where go :: (Int, Index) -> Maybe (Char, (Int, Index))
        go (0, rest) | rest == 0 = Nothing
                     | otherwise = error $ "fromIndex: Residual " ++ show rest ++ " from " ++ show (len, idx)
        go (l, x) = let (x', ch) = x `divMod` 4
                    in Just (case ch of
                              0 -> 'A'
                              1 -> 'C'
                              2 -> 'G'
                              3 -> 'T'
                              _ -> error $ "fromIndex: Bad result " ++ show (x', ch) ++ " from " ++ show x ++ " divMod 4"
                             , ((l - 1), x'))
                        

maxIndex :: Int -> Index
maxIndex len = (4 ^ len) - 1

lenIndex :: Int -> Int
lenIndex len = 4 ^ len

countSeq :: UM.IOVector Int -> BS.ByteString -> IO ()
countSeq ctvec sequ | UM.length ctvec >= lenIndex (BS.length sequ)
                      = let idx = toIndex sequ
                        in UM.read ctvec (unIndex idx) >>= UM.write ctvec (unIndex idx) . (succ $!)
                    | otherwise = fail $ "Sequence " ++ show sequ ++ " length " ++ show (BS.length sequ) ++ " needs "
                                  ++ (show . lenIndex . BS.length $ sequ) ++ " but count vector length " ++ show (UM.length ctvec)
