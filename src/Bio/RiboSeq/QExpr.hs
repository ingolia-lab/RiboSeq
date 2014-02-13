{-# LANGUAGE OverloadedStrings #-}
module Bio.RiboSeq.QExpr
       where

import Control.Applicative
import Data.List (transpose)
import Control.Monad
import qualified Data.ByteString.Char8 as BS

import qualified Data.Attoparsec as AP (parseOnly)
import qualified Data.Attoparsec.Char8 as AP
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import Statistics.Quantile

data QExpr = QExpr { name :: !BS.ByteString
                   , size :: !Int
                   , count :: !Int
                   } deriving (Show, Eq, Ord)
                              
dens :: QExpr -> Double
dens qe = (fromIntegral $ count qe) / (fromIntegral $ size qe)

toLine :: QExpr -> BS.ByteString
toLine qe = BS.intercalate (BS.singleton '\t') 
            [ name qe
            , BS.pack . show . size $ qe
            , BS.pack . show . count $ qe
            ]
            
line :: AP.Parser QExpr
line = (QExpr <$>
        (AP.takeWhile (/= '\t') <* AP.char '\t' AP.<?> "gene name") <*>
        (AP.signed AP.decimal <* AP.char '\t' AP.<?> "gene size") <*>
        (AP.decimal <* AP.endOfLine AP.<?> "gene count")) AP.<?> "gene qexpr"

read :: FilePath -> IO [QExpr]
read = BS.readFile >=> 
       either (ioError . userError) return .
       AP.parseOnly (AP.manyTill line AP.endOfInput)

data FeatureQE = FeatureQE { feature :: !BS.ByteString
                           , fsize :: !Int
                           , quant :: !(U.Vector Int)
                           } deriving (Show)

data QETable = QETable { samples :: !(V.Vector BS.ByteString)
                       , features :: !(V.Vector FeatureQE)
                       }

headerstr :: BS.ByteString
headerstr = "# gene\tsize"

readTable :: FilePath -> IO QETable
readTable = liftM parse . BS.readFile

parse :: BS.ByteString -> QETable
parse qtstr = QETable { samples = samps, features = feats }
  where (sampleLine:featureLines) = BS.lines qtstr
        samps = either error id . AP.parseOnly header $ sampleLine
        nsamps = V.length samps
        feats = V.fromList $! map feat featureLines
          where feat l = either badfeat id . AP.parseOnly (featureq nsamps) $ l
                  where badfeat err = error $ "Bad feature: " ++ err ++ ": " ++ show l
        
header :: AP.Parser (V.Vector BS.ByteString)
header = AP.string headerstr *> (V.fromList <$> AP.manyTill sample AP.endOfInput)
  where sample = BS.copy <$> (AP.char '\t' *> AP.takeTill AP.isSpace)

featureq :: Int -> AP.Parser FeatureQE
featureq nsamp = do f <- {-# SCC "f" #-} BS.copy <$> AP.takeTill (== '\t')
                    sz <- {-# SCC "s" #-} AP.char '\t' *> AP.decimal
                    qs <- parseValues nsamp U.empty
                    return $! FeatureQE f sz qs
  where parseValues 0 v0 = AP.endOfInput *> pure v0
        parseValues n v0 = (AP.char '\t' *> (U.snoc v0 <$> AP.decimal)) >>= parseValues (n - 1)

unparses :: QETable -> String -> String
unparses = unparsesWith (\x -> shows x)

unparsesWith :: (Int -> String -> String) -> QETable -> String -> String
unparsesWith f qt = unparsesHeader . unparsesFeatqs
  where unparsesHeader = (BS.unpack headerstr ++) . unparsesSamples . ('\n' :)
        unparsesSamples = V.foldl' (\i sz -> i . unparsesSample sz) id (samples qt)
        unparsesSample sz = ('\t' :) . (BS.unpack sz ++)
        unparsesFeatqs = V.foldl' (\i fq -> i . unparsesFeatq fq) id (features qt)
        unparsesFeatq fq = (BS.unpack (feature fq) ++) . ('\t' :) .
                           (shows $ fsize fq) .
                           U.foldl' unparseQuant id (quant fq) . ('\n' :)
        unparseQuant i q = i . ('\t' :) . f q

sampleQuants :: [[QExpr]] -> [(Double, Double)]
sampleQuants samples = zip normedquants factors
  where rawquants = map samplequant . transpose . filter (not . allZero) $ samples
        allZero = all ((== 0) . count)
        samplequant = continuousBy medianUnbiased 3 4 . U.fromList . map (fromIntegral . count)
        totalCounts = sum . map (fromIntegral . count) . concat $ samples
        totalQuants = sum rawquants
        normedquants = map (* (totalCounts / totalQuants)) rawquants
        factors = map (1e6 /) normedquants
