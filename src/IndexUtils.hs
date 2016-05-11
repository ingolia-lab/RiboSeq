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
import qualified Data.HashSet as HS
import Data.IORef
import Data.List
import Data.Maybe
import Data.Ord
import Numeric
import System.Directory
import System.Exit
import System.FilePath
import System.IO

import System.Console.CmdTheLine

main :: IO ()
main = runChoice ( indexUtils <$> baseConf, baseInfo )
       [ ( distances <$> distConf, distInfo )
       , ( generate <$> genConf, genInfo )
       , ( select <$> selConf, selInfo )
       ]
  where baseInfo = defTI { termName = "index-utils"
                         , version = "0.0"
                         , termDoc = "Utilities for index sequences"
                         }
        distInfo = defTI { termName = "distance"
                         , version = "0.0"
                         , termDoc = "Distance matrix of indexes"
                         }
        genInfo = defTI { termName = "generate"
                        , version = "0.0"
                        , termDoc = "Generate table of indexes"
                        }
        selInfo = defTI{ termName = "select"
                       , version = "0.0"
                       , termDoc = "Select non-conflicting indexes"
                       }

indexUtils :: BaseConf -> IO ()
indexUtils _ = return ()

data BaseConf = BaseConf

baseConf :: Term BaseConf
baseConf = pure BaseConf

distances :: DistConf -> IO ()
distances conf = BS.readFile (distInput conf) >>= BS.writeFile (distOutput conf) . distanceTable . parseIndex

distanceTable :: [(BS.ByteString, BS.ByteString)] -> BS.ByteString
distanceTable indexseqs = BS.unlines $ header : (map distanceLine indexseqs)
  where header = unfields $ (BS.pack "#") : (map fst indexseqs)
        distanceLine (name, sequ) = unfields $ name : distanceFields ++ worstFields
          where distances = map (\(n,s) -> (n, sequ `indexDist` s)) indexseqs
                distanceFields = map (BS.pack . show . snd) distances
                (worstName, worstDist) = minimumBy (comparing snd) . filter ((/= name) . fst) $ distances
                worstFields = [ worstName, BS.pack . show $ worstDist ]
        unfields = BS.intercalate (BS.singleton '\t')
  
data DistConf = DistConf { distInput :: FilePath
                         , distOutput :: FilePath
                         } deriving (Read, Show)

distConf :: Term DistConf
distConf = DistConf <$>
           input <*>
           output
  where input :: Term FilePath
        input = required $ opt Nothing $ (optInfo [ "i", "input" ])
                { optName = "INPUT", optDoc = "Input text file of index sequences" }
        output :: Term FilePath
        output = required $ opt Nothing $ (optInfo [ "o", "output" ])
                 { optName = "OUTPUT", optDoc = "Output table of index distances" }

-- Generate

generate :: GenConf -> IO ()
generate conf = BS.writeFile (genOutput conf) table
  where table = BS.unlines . zipWith indexLine [0..] $ indexSeqs
        indexSeqs = genNext (genLength conf) ""
        indexLine n sequ = BS.pack $ (show n) ++ "\t" ++ sequ

genNext :: Int -> String -> [String]
genNext 0 str = [ str ]
genNext l [] = concatMap (genNext (l - 1)) [ "A", "C", "G", "T" ]
genNext l str@[nt] = concatMap (genNext (l - 1) . (: str)) $ filter (/= nt) "ACGT"
genNext l str@(nt1:nt2:_) = concatMap (genNext (l - 1) . (: str)) .
                            filter (/= nt1) . filter (/= nt2) $ "ACGT"

data GenConf = GenConf { genLength :: Int
                       , genOutput :: FilePath
                       } deriving (Read, Show)

genConf :: Term GenConf
genConf = GenConf <$>
          len <*>
          output
  where len :: Term Int
        len  = required $ opt Nothing $ (optInfo [ "l", "length" ])
               { optName = "LENGTH", optDoc = "Length of index sequences" }
        output :: Term FilePath
        output = required $ opt Nothing $ (optInfo [ "o", "output" ])
                 { optName = "OUTPUT", optDoc = "Output table of index sequences" }

-- Select

select :: SelConf -> IO ()
select conf = do idxs <- liftM parseIndex $ BS.readFile (selInput conf)
                 exts <- case selExtant conf of
                   Nothing -> return []
                   Just extfile -> liftM parseIndex $ BS.readFile extfile
                 let !sel = selectify conf idxs exts
                 case selDetails conf of
                   Nothing -> return ()
                   Just detfile -> BS.writeFile detfile $ selectedDetails sel
                 BS.writeFile (selOutput conf) $ selectedTable sel                 

data Selection = Selection { sDist :: HM.HashMap (BS.ByteString, BS.ByteString) Int
                           , sSeqs :: HM.HashMap BS.ByteString BS.ByteString
                           , sPicked :: [(BS.ByteString, [BS.ByteString])]
                           , sActive :: HS.HashSet BS.ByteString
                           }

selectify :: SelConf -> [(BS.ByteString, BS.ByteString)] -> [(BS.ByteString, BS.ByteString)] -> Selection
selectify conf idxs exts = iterate $ newSelection conf idxs exts
  where iterate sel0 = maybe sel0 (iterate . selectIndex conf sel0) $ nextSelection conf sel0

newSelection :: SelConf -> [(BS.ByteString, BS.ByteString)] -> [(BS.ByteString, BS.ByteString)] -> Selection
newSelection conf idxs exts = Selection { sDist = dists, sSeqs = seqs, sPicked = [], sActive = HS.fromList $ map fst allIndexes }
  where allIndexes = idxs ++ exts
        dists = HM.fromList [ ((j, k), indexDist js ks) | (j, js) <- allIndexes, (k, ks) <- allIndexes ]
        seqs = HM.fromList allIndexes

nextSelectionDumb :: SelConf -> Selection -> Maybe BS.ByteString
nextSelectionDumb conf sel = case HS.toList (sActive sel) of
  [] -> Nothing
  (k0:rest) -> Just k0

nextSelection :: SelConf -> Selection -> Maybe BS.ByteString
nextSelection conf sel = case HS.toList (sActive sel) of
  [] -> Nothing
  l@(_:_) -> Just $! minimumBy (comparing nInactivated) l
  where nInactivated k = HS.size . HS.filter inactivated $ sActive sel
          where inactivated j = maybe True (< (selMinDist conf)) $
                                HM.lookup (j, k) (sDist sel)

selectIndex :: SelConf -> Selection -> BS.ByteString -> Selection
selectIndex conf sel0 k = sel0 { sPicked = (k, inactives) : (sPicked sel0)
                               , sActive = actives
                               }
  where actives = HS.filter stillActive $ sActive sel0
        inactives = filter (not . stillActive) . HS.toList . HS.delete k $ sActive sel0
        stillActive j = maybe False (>= (selMinDist conf)) $
                        HM.lookup (j, k) (sDist sel0)

selectedTable :: Selection -> BS.ByteString
selectedTable sel = BS.unlines . map pickedLine . reverse . sPicked $ sel
  where pickedLine (k, _) = BS.concat [ k, "\t", HM.lookupDefault "!!!" k $ sSeqs sel ]

selectedDetails :: Selection -> BS.ByteString
selectedDetails sel = BS.unlines . map pickedLine . reverse . sPicked $ sel
  where pickedLine (k, es) = BS.unwords $ (sequ k) : "excluded" : (map sequ es)
        sequ k = BS.concat [ k, "(", HM.lookupDefault "!!!" k $ sSeqs sel, ")" ]

data SelConf = SelConf { selMinDist :: Int
                       , selInput :: FilePath
                       , selExtant :: (Maybe FilePath)
                       , selOutput :: FilePath
                       , selDetails :: (Maybe FilePath)
                       } deriving (Read, Show)

selConf :: Term SelConf
selConf = SelConf <$> mindist <*> input <*> extant <*> output <*> details
  where mindist :: Term Int
        mindist = required $ opt Nothing $ (optInfo [ "m", "mindist" ])
                  { optName = "DIST", optDoc = "Minimum mismatches between index sequences" }
        input :: Term FilePath
        input = required $ opt Nothing $ (optInfo [ "i", "input" ])
                { optName = "INPUT", optDoc = "Input text file of index sequences" }
        extant :: Term (Maybe FilePath)
        extant = value $ opt Nothing $ (optInfo [ "e", "extant" ])
                 { optName = "EXTANT", optDoc = "Input text file of extant (required) index sequences" }
        output :: Term FilePath
        output = required $ opt Nothing $ (optInfo [ "o", "output" ])
                 { optName = "OUTPUT", optDoc = "Output table of index selections" }
        details :: Term (Maybe FilePath)
        details = value $ opt Nothing $ (optInfo [ "d", "details" ])
                 { optName = "DETAILS", optDoc = "Output details of index selection and exclusion" }
                  
--------

parseIndex :: BS.ByteString -> [(BS.ByteString, BS.ByteString)]
parseIndex = map parseIndexseq . filter notComment . BS.lines 
  where notComment = not . BS.isPrefixOf (BS.singleton '#')
        parseIndexseq l = case BS.split '\t' l of
          [ name, sequ ] -> (name, sequ)
          _ -> error $ "Malformed indexseq line " ++ show l
        
indexDist :: BS.ByteString -> BS.ByteString -> Int
indexDist s1 = sum . BS.zipWith (\c1 c2 -> if toUpper c1 == toUpper c2 then 0 else 1) s1
