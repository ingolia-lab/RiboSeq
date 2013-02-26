{-# LANGUAGE OverloadedStrings #-}

module Main
       where 

import Control.Applicative
import Control.Arrow
import Control.Exception
import Control.Monad.Reader
import qualified Data.ByteString.Char8 as BS
import Data.IORef
import Data.List
import Data.Maybe
import Numeric
import System.Console.GetOpt
import System.Environment
import System.FilePath
import System.IO

import qualified Data.Iteratee.IO as IterIO
import qualified Data.Iteratee.ListLike as IterLL

import qualified Bio.SamTools.BamIndex as BamIndex
import qualified Bio.SamTools.Bam as Bam
import qualified Bio.SamTools.Iteratee as BamIter
import Bio.SeqLoc.Bed
import Bio.SeqLoc.LocRepr
import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.Position as Pos
import Bio.SeqLoc.Strand
import Bio.SeqLoc.Transcript

import Bio.RiboSeq.CodonAssignment

verbose :: Bool
verbose = False

main :: IO ()
main = getArgs >>= handleOpt . getOpt RequireOrder optDescrs
    where handleOpt (_,    _,         errs@(_:_)) = usage (unlines errs)
          handleOpt (args, [bam], []) = either usage (doCount bam) $ argsToConf args
          handleOpt (_,    _,     []) = usage "Specify just one sorted, indexed BAM file"
          usage errs = do prog <- getProgName
                          let progline = prog ++ " [OPTIONS] <BAM>"
                          hPutStr stderr $ usageInfo progline optDescrs
                          hPutStrLn stderr errs

doCount :: FilePath -> Conf -> IO ()
doCount bam conf = (maybe doBamCount doBedCount $ confBed conf) bam conf

doBamCount :: FilePath -> Conf -> IO ()
doBamCount bam conf = do ct <- newCount conf
                         IterLL.run $ IterLL.joinIM $ BamIter.enumBam bam $ IterLL.mapM_ $ countBam ct
                         ctstr <- showCount ct
                         putStrLn $ intercalate "\t" [ bam, ctstr ]
                           
doBedCount :: FilePath -> FilePath -> Conf -> IO ()
doBedCount bed bam conf = do ct <- newCount conf
                             BamIndex.withIndex bam $ \hidx ->
                               let trxIter = bedTranscriptEnum . IterLL.mapM_ $ countTrxReads hidx ct 
                               in IterIO.fileDriver trxIter bed
                             ctstr <- showCount ct
                             putStrLn $ intercalate "\t" [ bam, bed, ctstr ]

countTrxReads :: BamIndex.IdxHandle -> Count -> Transcript -> IO ()
countTrxReads hidx ct trx = maybe noTarget withTarget $ Bam.lookupTarget header $ unSeqLabel trxref
  where header = BamIndex.idxHeader hidx
        (OnSeq trxref trxsploc) = location trx
        trxbnds = (fromIntegral *** fromIntegral) $ Loc.bounds trxsploc
        noTarget = when verbose $ hPutStrLn stderr . unwords $ [ "Skipping", show . unSeqLabel . geneId $ trx
                                                               , "on", show . unSeqLabel $ trxref ]
        withTarget tid = do when verbose $ hPutStrLn stderr . unwords $ [ "Running", show . unSeqLabel . geneId $ trx ]
                            IterLL.run targetIter
          where targetIter = IterLL.joinIM $ BamIter.enumIndexRegion hidx tid trxbnds $ IterLL.mapM_ countRead
        countRead b = case Bam.refSpLoc b >>= trxReadContig trx of
          Just _cloc -> countBam ct b
          _ -> return ()

data Count = Unscaled !(IORef Int)
           | Scaled !(IORef Double)
             
newCount :: Conf -> IO Count
newCount conf | confWeightNH conf = liftM Scaled $! newIORef 0.0
              | otherwise         = liftM Unscaled $! newIORef 0
                           
countBam :: Count -> Bam.Bam1 -> IO ()                                    
countBam (Unscaled zref) bam = do z0 <- readIORef zref
                                  writeIORef zref $! succ z0
countBam (Scaled xref) bam = maybe (return ()) count $! Bam.nHits bam
  where count nh = do x0 <- readIORef xref
                      writeIORef xref $! x0 + (recip $ fromIntegral nh)
                           
showCount :: Count -> IO String
showCount (Unscaled zref) = liftM show $! readIORef zref
showCount (Scaled xref) = do x <- readIORef xref
                             return $! showFFloat (Just 1) x ""

data Conf = Conf { confBed :: !(Maybe FilePath)
                 , confWeightNH :: !Bool
                 } deriving (Show)

data Arg = ArgBed { unArgBed :: !String }
         | ArgASite { unArgASite :: !String }
         | ArgWeightNH
         deriving (Show, Read, Eq, Ord)

argBed :: Arg -> Maybe String
argBed (ArgBed bed) = Just bed
argBed _ = Nothing

argASite :: Arg -> Maybe String
argASite (ArgASite aSite) = Just aSite
argASite _ = Nothing

optDescrs :: [OptDescr Arg]
optDescrs = [ Option ['b'] ["bed"]        (ReqArg ArgBed "BED")          "Bed filename"
            , Option ['w'] ["weight-nh"]  (NoArg ArgWeightNH)            "Weight reads by NH"
            ]

argsToConf :: [Arg] -> Either String Conf
argsToConf = runReaderT conf
    where conf = Conf <$> 
                 findBed <*>
                 ReaderT (return . elem ArgWeightNH)
          findBed = ReaderT $ return . listToMaybe . mapMaybe argBed

