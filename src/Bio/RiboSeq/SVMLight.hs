{-# LANGUAGE BangPatterns #-}
module Bio.RiboSeq.SVMLight
    where

import Control.Monad
import qualified Data.ByteString.Char8 as BS
import Data.Char (isSpace)
import Data.List
import Data.Monoid
import Numeric
import System.Exit
import System.FilePath
import System.IO (IOMode(..), withFile)
import System.Process

import qualified Data.Vector.Unboxed as U

verbose :: Bool
verbose = True

myRawSystem :: String -> [String] -> IO ExitCode
myRawSystem execu arglist = do
  when verbose $ putStrLn (unwords $ execu : arglist)
  rawSystem execu arglist

targetYes, targetNo, targetQuery :: BS.ByteString
targetYes = BS.pack "+1"
targetNo = BS.pack "-1"
targetQuery = BS.pack "0"

exampleLine :: BS.ByteString -> [(Int, Double)] -> Maybe BS.ByteString -> BS.ByteString
exampleLine target features minfo = BS.unwords $ target : (BS.pack featuresStr) : infobs
    where infobs = maybe [] (\info -> [ BS.singleton '#', info ]) minfo
          featuresStr = foldr ($) "" . intersperse (' ' :) . map showsFeature $ features
          showsFeature (fidx, fval) = shows fidx . (':' :) . showEFloat (Just 4) fval

data SVMMode = SVMClassification | SVMRegression | SVMPreference deriving (Show, Eq, Ord, Enum)

modeOpts :: SVMMode -> [String]
modeOpts SVMClassification = [ "-z", "c" ]
modeOpts SVMRegression     = [ "-z", "r" ]
modeOpts SVMPreference     = [ "-z", "p" ]

data SVMKernel = SVMLinear
               | SVMPolynomial { degree :: Maybe Int, scaling :: Maybe Double, baseline :: Maybe Double }
               | SVMRadial { gamma :: Maybe Double }
               | SVMSigmoid { scaling :: Maybe Double, baseline :: Maybe Double }
                 deriving (Show)

kernelOpts :: SVMKernel -> [String]
kernelOpts SVMLinear = [ "-t", "0" ]
kernelOpts (SVMPolynomial md ms mc) = concat [ [ "-t", "1" ]
                                             , maybe [] (\d -> [ "-d", show d ]) md
                                             , maybe [] (\s -> [ "-s", show s ]) ms
                                             , maybe [] (\c -> [ "-c", show c ]) mc
                                             ]
kernelOpts (SVMRadial mg) = [ "-t", "2" ] ++ maybe [] (\g -> [ "-g", show g ]) mg
kernelOpts (SVMSigmoid ms mc) = concat [ [ "-t", "3" ]
                                       , maybe [] (\s -> [ "-s", show s ]) ms
                                       , maybe [] (\c -> [ "-c", show c ]) mc
                                       ]

data RunSVMLearn = RunSVMLearn { learnBinary :: FilePath
                               , examplesIn :: FilePath
                               , modelOut :: FilePath
                               , mode :: Maybe SVMMode          -- -z
                               , learnVerbose :: Maybe Int      -- -v
                               , softMargin :: Maybe Double     -- -c
                               , width :: Maybe Double          -- -w
                               , positiveWeight :: Maybe Double -- -j
                               , biasedPlane :: Maybe Bool      -- -b
                               , removeRetrain :: Maybe Bool    -- -i
                               , positiveTrans :: Maybe Double  -- -p
                               , kernel :: Maybe SVMKernel      -- -t and -d, -g, -s, -r, -u
                               , maxQPSize :: Maybe Int         -- -q
                               , maxNewVar :: Maybe Int         -- -n
                               , kernelCache :: Maybe Int       -- -m
                               , termEps :: Maybe Double        -- -e
                               , shrinkIter :: Maybe Int        -- -h
                               , checkFinalOpt :: Maybe Bool    -- -f
                               , maxIter :: Maybe Int           -- -#
                               } deriving (Show)

defaultLearn :: RunSVMLearn
defaultLearn = RunSVMLearn { learnBinary = error "SVMLight: no svm_learn binary specified"
                           , examplesIn = error "SVMLight: no examples input file specified"
                           , modelOut = error "SVMLight: no model output file specified"
                           , mode = Nothing, learnVerbose = Nothing, softMargin = Nothing, width = Nothing
                           , positiveWeight = Nothing, biasedPlane = Nothing, removeRetrain = Nothing
                           , positiveTrans = Nothing
                           , kernel = Nothing
                           , maxQPSize = Nothing, maxNewVar = Nothing, kernelCache = Nothing, termEps = Nothing
                           , shrinkIter = Nothing, checkFinalOpt = Nothing, maxIter = Nothing
                           }

learnOpts :: RunSVMLearn -> [String]
learnOpts svmLearn = concat [ zOpts, vOpts, cOpts, wOpts, jOpts, bOpts, iOpts
                            , pOpts
                            , maybe [] kernelOpts . kernel $ svmLearn 
                            , qOpts, nOpts, mOpts, eOpts, hOpts, fOpts, poundOpts
                            ]
    where zOpts = maybe [] modeOpts . mode $ svmLearn
          vOpts = maybe [] (\v -> [ "-v", show v ]) . learnVerbose $ svmLearn
          cOpts = maybe [] (\c -> [ "-c", show c ]) . softMargin $ svmLearn
          wOpts = maybe [] (\w -> [ "-w", show w ]) . width $ svmLearn
          jOpts = maybe [] (\j -> [ "-j", show j ]) . positiveWeight $ svmLearn
          bOpts = maybe [] (\b -> [ "-b", if b then "1" else "0" ]) . biasedPlane $ svmLearn
          iOpts = maybe [] (\i -> [ "-i", if i then "1" else "0" ]) . removeRetrain $ svmLearn
          pOpts = maybe [] (\p -> [ "-p", show p ]) . positiveTrans $ svmLearn
          qOpts = maybe [] (\q -> [ "-q", show q ]) . maxQPSize $ svmLearn
          nOpts = maybe [] (\n -> [ "-n", show n ]) . maxNewVar $ svmLearn
          mOpts = maybe [] (\m -> [ "-m", show m ]) . kernelCache $ svmLearn
          eOpts = maybe [] (\e -> [ "-e", show e ]) . termEps $ svmLearn
          hOpts = maybe [] (\h -> [ "-h", show h ]) . shrinkIter $ svmLearn
          fOpts = maybe [] (\f -> [ "-f", if f then "1" else "0" ]) . checkFinalOpt $ svmLearn
          poundOpts = maybe [] (\p -> [ "-#", show p ]) . maxIter $ svmLearn

runSVMLearn :: RunSVMLearn -> IO ExitCode
runSVMLearn svmLearn = myRawSystem (learnBinary svmLearn) $ learnOpts svmLearn ++ [ examplesIn svmLearn, modelOut svmLearn ]

data RunSVMClassify = RunSVMClassify { classifyBinary :: FilePath
                                     , queriesIn :: FilePath
                                     , modelIn :: FilePath
                                     , decisionsOut :: FilePath
                                     , classifyVerbose :: Maybe Int -- -v
                                     } deriving (Show)

defaultClassify :: RunSVMClassify
defaultClassify = RunSVMClassify { classifyBinary = error "SVMLight: no svm_classify binary specified"
                                 , queriesIn = error "SVMLight: no examples input file specified"
                                 , modelIn = error "SVMLight: no model input file specified"
                                 , decisionsOut = error "SVMLight: No decisions output file specified"
                                 , classifyVerbose = Nothing
                                 }

classifyOpts :: RunSVMClassify -> [String]
classifyOpts svmLearn = concat [ vOpts ]
    where vOpts = maybe [] (\v -> [ "-v", show v ]) . classifyVerbose $ svmLearn

runSVMClassify :: RunSVMClassify -> IO ExitCode
runSVMClassify svmClassify = myRawSystem (classifyBinary svmClassify) $
                             classifyOpts svmClassify ++ [ queriesIn svmClassify, modelIn svmClassify, decisionsOut svmClassify ]

-- runSVMClassifyWithScore :: RunSVMClassify -> IO ExitCode
-- runSVMClassifyWithScore classify = bracket mkScoreTmp rmScoreTmp $ \scorename -> 
--                                    let classifyTmp = classify { decisionsOut = scorename }
--                                    in runSVMClassify classifyTmp >>= \ec -> do
--                                        when (ec == ExitSuccess) $ 
--                                             pasteScore (queriesIn classify) scorename (decisionsOut classify)
--                                        return ec
--     where mkScoreTmp = mkstemp (decisionsOut classify ++ "XXXXXX") >>= \(name, h) ->
--                        hClose h >> return name
--           rmScoreTmp scorename = doesFileExist scorename >>= flip when (removeFile scorename)

pasteScore :: FilePath -> FilePath -> FilePath -> IO ()
pasteScore queryFilename scoreFilename pastedFilename
    = do query <- liftM BS.lines . BS.readFile $ queryFilename
         score <- liftM BS.lines . BS.readFile $ scoreFilename 
         withFile pastedFilename WriteMode $ \hout ->
             mapM_ (BS.hPutStrLn hout) $ map pasteLine $ mergeScore query score
  where base = BS.dropWhile isSpace . BS.drop 1
        unscored q0 = BS.concat [ base q0, BS.pack "\tn/a" ]
        scored q0 sc0 = BS.concat [ base . BS.dropWhile (/= '#') $ q0, BS.singleton '\t', sc0 ]
        pasteLine (q0, sc0) | BS.null sc0 = unscored q0
                            | otherwise = scored q0 sc0

mergeScore :: [BS.ByteString] -> [BS.ByteString] -> [(BS.ByteString, BS.ByteString)]
mergeScore query score = unfoldr merge (query, score)
  where merge ([], []) = Nothing
        merge ([], rest) = error $ "pasteScore: remaining scores: " ++ show rest
        merge ((q0:qrest), [])
          | isComment q0 = Just ((q0, BS.empty), (qrest, []))
          | otherwise = error $ "Ran out of scores for " ++ show q0
        merge ((q0:qrest), scs@(sc0:screst))
          | isComment q0 = Just ((q0, BS.empty), (qrest, scs))
          | otherwise    = Just ((q0, sc0),      (qrest, screst))
        isComment = maybe blankQuery ((== '#') . fst) . BS.uncons
        blankQuery = error "Blank example line"

data TestData a = TestData {  truePos, falsePos, trueNeg, falseNeg :: !a } deriving (Show)
type TestCount = TestData Int
type TestScore = TestData (U.Vector Double)

instance Monoid a => Monoid (TestData a) where
  mempty = TestData mempty mempty mempty mempty
  (TestData tp1 fp1 tn1 fn1) `mappend` (TestData tp2 fp2 tn2 fn2) 
    = TestData (tp1 `mappend` tp2) (fp1 `mappend` fp2) (tn1 `mappend` tn2) (fn1 `mappend` fn2)
  mconcat tds = TestData { truePos = mconcat $ map truePos tds
                         , falsePos = mconcat $ map falsePos tds
                         , trueNeg = mconcat $ map trueNeg tds
                         , falseNeg = mconcat $ map falseNeg tds
                         }
  
negatives :: TestCount -> Int
negatives tc = trueNeg tc + falsePos tc

positives :: TestCount -> Int
positives tc = truePos tc + falseNeg tc

calledPositives :: TestCount -> Int
calledPositives tc = truePos tc + falsePos tc

calledNegatives :: TestCount -> Int
calledNegatives tc = trueNeg tc + falseNeg tc

specificity :: TestCount -> Double
specificity tc = (fromIntegral . trueNeg $ tc) / (fromIntegral . negatives $ tc)

sensitivity :: TestCount -> Double
sensitivity tc = (fromIntegral . truePos $ tc) / (fromIntegral . positives $ tc)

displayTestCount :: TestCount -> String
displayTestCount tc = intercalate "\t" [ show . truePos $ tc, show . falsePos $ tc
                                       , show . trueNeg $ tc, show . falseNeg $ tc
                                       , showf . specificity $ tc
                                       , showf . sensitivity $ tc
                                       ]
    where showf = flip (showFFloat (Just 3)) ""

checkTest :: FilePath -> FilePath -> IO TestCount
checkTest examples score = do ex <- liftM BS.lines . BS.readFile $ examples
                              sc <- liftM BS.lines . BS.readFile $ score
                              let ls = filter scored $ mergeScore ex sc
                              return $! foldl' addScore (TestData 0 0 0 0) ls
    where addScore cts0 (exline, scline) = let exval = readSignedDouble . BS.takeWhile (not . isSpace) $ exline
                                               scval = readSignedDouble . BS.takeWhile (not . isSpace) $ scline
                                           in case (exval > 0, scval > 0) of
                                                (True,  True ) -> cts0 { truePos = truePos cts0 + 1 }
                                                (False, True ) -> cts0 { falsePos = falsePos cts0 + 1 }
                                                (False, False) -> cts0 { trueNeg = trueNeg cts0 + 1 }
                                                (True,  False) -> cts0 { falseNeg = falseNeg cts0 + 1 }
          readSignedDouble :: BS.ByteString -> Double
          readSignedDouble xstr = case BS.uncons xstr of
                                    Just ('+', pxstr) -> readDouble pxstr
                                    _ -> readDouble xstr
          readDouble xstr = case reads . BS.unpack $ xstr of
                              [(x, "")] -> x
                              _ -> error $ "readDouble: malformed Double " ++ show xstr
          scored (_q0, sc0) = not . BS.null $ sc0

partitionTest :: FilePath -> FilePath -> IO TestScore
partitionTest examples score = do ex <- liftM BS.lines . BS.readFile $ examples
                                  sc <- liftM BS.lines . BS.readFile $ score
                                  let ls = filter scored $ mergeScore ex sc
                                  return $! foldl' addScore (TestData U.empty U.empty U.empty U.empty) ls
    where addScore tsc0 (exline, scline) 
            = let exval = readSignedDouble . BS.takeWhile (not . isSpace) $ exline
                  scval = readSignedDouble . BS.takeWhile (not . isSpace) $ scline
              in case (exval > 0, scval > 0) of
                (True,  True ) -> tsc0 { truePos = truePos tsc0 `U.snoc` scval }
                (False, True ) -> tsc0 { falsePos = falsePos tsc0 `U.snoc` scval }
                (False, False) -> tsc0 { trueNeg = trueNeg tsc0 `U.snoc` scval }
                (True,  False) -> tsc0 { falseNeg = falseNeg tsc0 `U.snoc` scval }
          readSignedDouble :: BS.ByteString -> Double
          readSignedDouble xstr = case BS.uncons xstr of
                                    Just ('+', pxstr) -> readDouble pxstr
                                    _ -> readDouble xstr
          readDouble xstr = case reads . BS.unpack $ xstr of
                              [(x, "")] -> x
                              _ -> error $ "readDouble: malformed Double " ++ show xstr
          scored (_q0, sc0) = not . BS.null $ sc0

countPartition :: TestScore -> TestCount
countPartition tsc0 = TestData { truePos = U.length . truePos $ tsc0
                               , falsePos = U.length . falsePos $ tsc0
                               , trueNeg = U.length . trueNeg $ tsc0
                               , falseNeg = U.length . falseNeg $ tsc0
                               }