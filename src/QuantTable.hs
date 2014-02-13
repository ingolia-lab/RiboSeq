{-# LANGUAGE OverloadedStrings #-}

module Main
       where 

import Control.Applicative
import Control.Monad.Reader
import qualified Data.ByteString.Char8 as BS
import Data.List
import Data.Maybe
import qualified Data.Set as S
import Numeric
import System.Console.GetOpt
import System.Environment
import System.FilePath
import System.IO

import qualified Data.Attoparsec as AP (parseOnly)
import qualified Data.Attoparsec.Char8 as AP
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector as V
import Statistics.Quantile

import qualified Bio.RiboSeq.QExpr as QE

main :: IO ()
main = getArgs >>= handleOpt . getOpt RequireOrder optDescrs 
  where handleOpt (_,    _,         errs@(_:_)) = usage (unlines errs)
        handleOpt (args, [], []) = either usage doQuantTable $ argsToConf args
        handleOpt (_,    _,  []) = usage "Unexpected non-option arguments"
        usage errs = do prog <- getProgName
                        let progline = prog ++ " [OPTIONS]"
                        hPutStr stderr $ usageInfo progline optDescrs
                        hPutStrLn stderr errs

doQuantTable :: Conf -> IO ()
doQuantTable conf = do samples <- liftM transpose $! mapM (QE.read . snd) $ confQuantSample conf
                       wantedGenes <- readIsInIdList conf
                       let wantedSamples = mapMaybe (wantedFeature wantedGenes) samples
                       when (confQuantNorm conf) $ writeQuantFactors conf wantedSamples
                       writeTable conf wantedSamples
                        
readIsInIdList :: Conf -> IO (BS.ByteString -> Bool)
readIsInIdList = maybe (return $ const True) (liftM isIn . BS.readFile) . confIdList
  where isIn idfile = flip S.member (S.fromList $ BS.lines idfile)

wantedFeature :: (BS.ByteString -> Bool) -> [QE.QExpr] -> Maybe [QE.QExpr]
wantedFeature _ [] = error "No samples"
wantedFeature wantedGene ss@(s0:_) | any mismatchKey ss = error . unwords $ "Mismatch:" : map show ss
                                   | wantedGene (QE.name s0) = Just ss
                                   | otherwise = Nothing
                                     where mismatchKey s' = (QE.name s' /= QE.name s0) ||
                                                            (QE.size s' /= QE.size s0)
                      
writeTable :: Conf -> [[QE.QExpr]] -> IO ()
writeTable conf genes = writeFile (confOutput conf) $ QE.unparses qt ""
  where qt = QE.QETable { QE.samples = V.fromList $ map (BS.pack . fst) $ confQuantSample conf
                        , QE.features = V.fromList $ map qtfeat genes
                        }
        qtfeat [] = error "No samples"
        qtfeat ss@(s0:_) = QE.FeatureQE { QE.feature = QE.name s0
                                        , QE.fsize = QE.size s0
                                        , QE.quant = U.fromList . map QE.count $ ss
                                        }

writeQuantFactors :: Conf -> [[QE.QExpr]] -> IO ()
writeQuantFactors conf wantedSamples = withFile (confQuantOut conf) WriteMode $ \hout ->
  let pquant name (q, f) = hPutStrLn hout qline
        where qline = intercalate "\t" [ name, showf 0 q, showf 1 (q / 1e6), showf 6 f ]
              showf sf x = showFFloat (Just sf) x ""
      qs = QE.sampleQuants wantedSamples
      names = map fst $ confQuantSample conf
  in sequence_ $! zipWith pquant names qs

data Conf = Conf { confOutput :: !(FilePath) 
                 , confQuantSample :: ![(String, FilePath)]
                 , confIdList :: !(Maybe FilePath)
                 , confQuantNorm :: !Bool
                 } deriving (Show)

confQuantOut :: Conf -> FilePath
confQuantOut conf = let (base, ext) = splitExtension . confOutput $ conf
                    in (base ++ "_qnorm") `addExtension` ext

data Arg = ArgOutput { unArgOutput :: !String }
         | ArgQuantSample { unArgQuantSample :: !String }
         | ArgIdList { unArgIdList :: !String }
         | ArgQuantNorm
         | ArgQuantDens
         deriving (Show, Read, Eq, Ord)

argOutput :: Arg -> Maybe String
argOutput (ArgOutput del) = Just del
argOutput _ = Nothing

argQuantSample :: Arg -> Maybe String
argQuantSample (ArgQuantSample quantCf) = Just quantCf
argQuantSample _ = Nothing

argIdList :: Arg -> Maybe String
argIdList (ArgIdList idList) = Just idList
argIdList _ = Nothing

optDescrs :: [OptDescr Arg]
optDescrs = [ Option ['o'] ["output"]    (ReqArg ArgOutput "OUTPUT")          "Output filename"
            , Option ['s'] ["sample"]    (ReqArg ArgQuantSample "NAME,QEXPR") "Comparison sample name and qexpr file"
            , Option ['l'] ["id-list"]   (ReqArg ArgIdList "FILENAME")        "Filename for list of gene IDs"
            , Option ['q'] ["quantnorm"] (NoArg ArgQuantNorm)                 "Quantile normalization"
            ]

argsToConf :: [Arg] -> Either String Conf
argsToConf = runReaderT conf
    where conf = Conf <$> 
                 findOutput <*>
                 findQuantSample <*>
                 findIdList <*>
                 (ReaderT $ return . elem ArgQuantNorm)
          findOutput = ReaderT $ maybe (Left "No out base") return  . listToMaybe . mapMaybe argOutput
          findQuantSample = ReaderT $ mapM parseQuant . mapMaybe argQuantSample
          findIdList = ReaderT $ return . listToMaybe . mapMaybe argIdList

parseQuant :: String -> Either String (String, FilePath)
parseQuant = AP.parseOnly quant . BS.pack
  where quant = ((,) <$> 
                 (BS.unpack <$> AP.takeWhile (/= ',') <* AP.char ',') <*> 
                 (BS.unpack <$> AP.takeByteString)) AP.<?> "quant sample"

unfields :: [String] -> String
unfields = intercalate "\t"