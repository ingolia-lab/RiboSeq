module Main
       where

import Control.Applicative
import Control.Exception
import Control.Monad
import Control.Monad.Reader
import qualified Data.ByteString.Char8 as BS
import qualified Data.HashMap.Strict as M
import Data.Maybe
import System.Console.GetOpt
import System.Environment
import System.IO

import qualified Data.Attoparsec as AP
import qualified Data.Attoparsec.Char8 as AP
import qualified Data.Iteratee as Iter

import qualified Bio.SamTools.Bam as Bam
import qualified Bio.SamTools.Iteratee as BamIter

import Bio.RiboSeq.RrnaStats
import Bio.RiboSeq.RrnaFrags

main :: IO ()
main = getArgs >>= handleOpt . getOpt RequireOrder optDescrs
    where handleOpt (_,    _,         errs@(_:_)) = usage (unlines errs)
          handleOpt (args, [bam], []) = either usage (doRrnaStats bam) $ argsToConf args
          handleOpt (_,    _,     []) = usage "Specify just one BAM file"
          usage errs = do prog <- getProgName
                          let progline = prog ++ " [OPTIONS] <BAM>"
                          hPutStr stderr $ usageInfo progline optDescrs
                          hPutStrLn stderr errs


doRrnaStats :: FilePath -> Conf -> IO ()
doRrnaStats bam conf = bracket (openInFile bam) Bam.closeInHandle $ \hin -> do
  statsio <- mkRrnaStatsIO (confMaxReadLen conf) (Bam.targetSeqList . Bam.inHeader $ hin)
  verbose conf $ "Read target sequences " ++ (show . M.keys . unRrnaStatsIO $ statsio)
  Iter.run =<< BamIter.enumInHandle hin (Iter.mapM_ $ countAlign statsio)
  stats <- freezeRrnaStatsIO statsio
  writeFile (confOutput conf ++ "_rrna_count.txt") $ countTable stats
  writeFile (confOutput conf ++ "_rrna_pos.txt") $ posTable stats
  writeFile (confOutput conf ++ "_rrna_poslen.txt") $ posLenTable (confLenRange conf) stats
  writeFile (confOutput conf ++ "_rrna_frags.txt") $ fragmentTable stats
  where openInFile | confTamInput conf = Bam.openTamInFile
                   | otherwise = Bam.openBamInFile

data Conf = Conf { confOutput :: !(FilePath) 
                 , confTamInput :: !Bool
                 , confMaxReadLen :: !Int
                 , confLenRange :: !(Int, Int)
                 , confVerbose :: !Bool
                 } deriving (Show)

defaultMaxReadLen :: Int
defaultMaxReadLen = 100

defaultLenRange :: (Int, Int)
defaultLenRange = (25, 34)

verbose :: Conf -> String -> IO ()
verbose conf = when (confVerbose conf) . hPutStrLn stderr

data Arg = ArgOutput { unArgOutput :: !String }
         | ArgIsTam
         | ArgMaxReadLen { unMaxReadLen :: !String }
         | ArgLenRange { unLenRange :: !String }
         | ArgVerbose
         deriving (Show, Read, Eq, Ord)

argOutput :: Arg -> Maybe String
argOutput (ArgOutput del) = Just del
argOutput _ = Nothing

argMaxReadLen :: Arg -> Maybe String
argMaxReadLen (ArgMaxReadLen maxReadLen) = Just maxReadLen
argMaxReadLen _ = Nothing

argLenRange :: Arg -> Maybe String
argLenRange (ArgLenRange lenRange) = Just lenRange
argLenRange _ = Nothing

optDescrs :: [OptDescr Arg]
optDescrs = [ Option ['o'] ["output"]   (ReqArg ArgOutput "OUTFILE")   "Output filename"
            , Option ['t'] ["tam"]      (NoArg ArgIsTam)               "Text-format TAM input"
            , Option ['m'] ["maxread"]  (ReqArg ArgMaxReadLen "LEN")   ("Maximum read length [" ++ show defaultMaxReadLen ++ "]")
            , Option ['l'] ["lenrange"] (ReqArg ArgLenRange "MIN,MAX") ("Length range for output tables [" ++ show defaultLenRange ++ "]")
            , Option ['v'] ["verbose"]  (NoArg ArgVerbose)             "Verbose"
            ]

argsToConf :: [Arg] -> Either String Conf
argsToConf = runReaderT conf
    where conf = Conf <$> 
                 findOutput <*>
                 findIsTam <*>
                 findMaxRead <*>
                 findLenRange <*>
                 findVerbose
          findOutput = ReaderT $ maybe (Left "No out base") return  . listToMaybe . mapMaybe argOutput
          findIsTam = ReaderT $ return . elem ArgIsTam
          findMaxRead = ReaderT $ maybe (return defaultMaxReadLen) parseInt . listToMaybe . mapMaybe argMaxReadLen
          findLenRange = ReaderT $ maybe (return defaultLenRange) parseLenRange . listToMaybe . mapMaybe argLenRange
          findVerbose = ReaderT $ return . elem ArgVerbose
          parseInt = AP.parseOnly (AP.decimal <* AP.endOfInput) . BS.pack
          parseLenRange = AP.parseOnly apLenRange . BS.pack
            where apLenRange = (,) <$> AP.decimal <*> (AP.char ',' *> AP.decimal <* AP.endOfInput)
 