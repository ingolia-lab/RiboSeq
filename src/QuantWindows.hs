{-# LANGUAGE BangPatterns #-}

module Main
       where

import Control.Applicative
import Control.Arrow
import Control.Exception
import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.Reader
import qualified Data.ByteString.Char8 as BS
import Data.IORef
import Data.Int
import Data.List
import Data.Maybe
import Numeric
import System.Console.GetOpt
import System.Environment
import System.FilePath
import System.IO

import qualified Data.Attoparsec.ByteString.Char8 as AP
import qualified Data.Iteratee.IO as IterIO
import qualified Data.Iteratee.ListLike as IterLL

import qualified Bio.SamTools.Bam as Bam
import qualified Bio.SamTools.BamIndex as BamIndex
import qualified Bio.SamTools.Iteratee as BamIter

import Bio.SeqLoc.Bed
import Bio.SeqLoc.LocRepr
import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.Position as Pos
import Bio.SeqLoc.Strand
import Bio.SeqLoc.Transcript

import Bio.RiboSeq.CodonAssignment

doCount :: FilePath -> FilePath -> Conf -> IO ()
doCount bamfile outfile conf = BamIndex.withIndex bamfile $ \hidx -> 
  withFile outfile WriteMode $ \hout ->  
  (maybe doBamCount doBedCount $ confBed conf) conf hidx hout

doBamCount :: Conf -> BamIndex.IdxHandle -> Handle -> IO ()
doBamCount conf hidx hout =
  forM_ [0..(Bam.nTargets (BamIndex.idxHeader hidx) - 1)] $ \targetid ->
  let seqlen = Bam.targetSeqLen (BamIndex.idxHeader hidx) targetid
      winlen = fromIntegral $ confWindow conf
      winstep = fromIntegral $ confWindowStep conf
      maxwin = 1 + ((seqlen - winlen) `div` winstep)
  in do when (confVerbose conf) $ do
          BS.hPutStr stderr $ Bam.targetSeqName (BamIndex.idxHeader hidx) targetid
          hFlush stderr
        forM_ [0..maxwin] $ \win ->
          let winstart = (winstep * win)
              winend = min (seqlen - 1) (winstart + winlen - 1)
          in targetWindow conf hidx targetid (winstart, winend) >>= hPutStr hout
        when (confVerbose conf) $ hPutStrLn stderr $ "..." ++ show maxwin
     
targetWindow :: Conf -> BamIndex.IdxHandle -> Int -> (Int64, Int64) -> IO String
targetWindow conf hidx targetid (wstart, wend) = query >>= flip BamIter.enumQuery window >>= IterLL.run
  where query = BamIndex.query hidx targetid (wstart, wend)
        seqname = Bam.targetSeqName (BamIndex.idxHeader hidx) targetid
        inwin b = maybe False (\o -> o >= wstart && o <= wend) $ Bam.position b
        window = do counts@(raw, _uniq, _perf, _scaled) <- IterLL.joinI $ IterLL.filter inwin
                                                           $ IterLL.zip4 iterRaw iterUniq iterPerf iterEffec  
                    if (raw >= confWindowMin conf) 
                      then return $ unlines [ intercalate "\t" $ [ concat [ BS.unpack seqname, ":", show wstart, "-", show wend ]
                                                                 , BS.unpack seqname
                                                                 , show wstart
                                                                 ] ++ countFields counts
                                            ]
                      else return ""

iterRaw, iterUniq, iterPerf :: (Functor m, Monad m) => IterLL.Iteratee [Bam.Bam1] m Int
iterRaw = IterLL.length

iterUniq = IterLL.joinI $ IterLL.filter isUniq IterLL.length
  where isUniq = maybe False (== 1) . Bam.nHits
        
iterPerf = IterLL.joinI $ IterLL.filter isPerf IterLL.length
  where isPerf = maybe False (== 0) . Bam.nMismatch

iterEffec :: (Monad m) => IterLL.Iteratee [Bam.Bam1] m Double
iterEffec = IterLL.joinI $ IterLL.mapStream effec IterLL.sum
  where effec = maybe 0 (recip . fromIntegral) . Bam.nHits

doBedCount :: FilePath -> Conf -> BamIndex.IdxHandle -> Handle -> IO ()
doBedCount bedfile conf hidx hout = IterIO.fileDriver trxIter bedfile
  where trxIter = bedTranscriptEnum . IterLL.mapM_ $ \trx ->
          countTrx conf hidx hout trx 

countTrx :: Conf -> BamIndex.IdxHandle -> Handle -> Transcript -> IO ()
countTrx conf hidx hout trx = maybe noTarget withTarget $ Bam.lookupTarget header $ unSeqLabel . onSeqLabel . location $ trx
  where header = BamIndex.idxHeader hidx
        noTarget = do when (confVerbose conf) $ hPutStrLn stderr . unwords $ 
                        [ "Skipping", show . unSeqLabel . geneId $ trx
                        , "on", show . unSeqLabel . onSeqLabel . location $ trx ]
        maxwin = Loc.length trxsploc `div` confWindow conf
        winpos w = w * (confWindow conf)
        (OnSeq trxref trxsploc) = location trx
        trxbnds = (fromIntegral *** fromIntegral) $ Loc.bounds trxsploc
        withTarget tid = do when (confVerbose conf) $ do
                              hPutStr stderr . unwords $ [ "Running", show . unSeqLabel . geneId $ trx ]
                              hFlush stderr
                            forM_ [0..maxwin] $ \win ->
                              (IterLL.run $ IterLL.joinIM $ BamIter.enumIndexRegion hidx tid trxbnds $
                               countTrxWindow conf trx (winpos win, winpos (win + 1) - 1)) >>=
                              hPutStr hout
                            when (confVerbose conf) $ hPutStrLn stderr $ "..." ++ show (fromIntegral maxwin :: Int)
                              
countTrxWindow :: (Monad m, MonadIO m, Functor m) => Conf -> Transcript -> (Pos.Offset, Pos.Offset) -> IterLL.Iteratee [Bam.Bam1] m String
countTrxWindow conf trx (offmin, offmax) 
  = let readStart b = liftM (Pos.offset . Loc.startPos) $ 
                      Bam.refSpLoc b >>= trxReadContig trx
        inWindow = maybe False (\o -> o >= offmin && o <= offmax) . readStart
        trxname = BS.unpack . unSeqLabel . geneId $ trx
    in do counts@(raw, _uniq, _perf, _scaled) <- IterLL.joinI $ IterLL.filter inWindow $ 
                                                 IterLL.zip4 iterRaw iterUniq iterPerf iterEffec  
          if (raw >= confWindowMin conf) 
             then return $ unlines [ intercalate "\t" $ [ concat [ trxname, ":"
                                                                 , show (fromIntegral offmin :: Int), "-"
                                                                 , show (fromIntegral offmax :: Int) 
                                                                 ]
                                                        , trxname
                                                        , show (fromIntegral offmin :: Int)
                                                        ] ++ countFields counts
                                   ]
             else return ""

countFields :: (Int, Int, Int, Double) -> [String]
countFields (raw, uniq, perf, scaled) = [ showFFloat (Just 1) scaled ""                                          
                                        , show raw
                                        , showFFloat (Just 3) (scaled / fromIntegral raw) ""
                                        , showFFloat (Just 3) (fromIntegral uniq / fromIntegral raw) ""
                                        , showFFloat (Just 3) (fromIntegral perf / fromIntegral raw) ""
                                        ]

main :: IO ()
main = getArgs >>= handleOpt . getOpt RequireOrder optDescrs
    where handleOpt (_,    _,          errs@(_:_)) = usage (unlines errs)
          handleOpt (args, [bam, out], []) = either usage (doCount bam out) $ argsToConf args
          handleOpt (_,    _,          []) = usage "Specify just one sorted, indexed BAM file and one output file"
          usage errs = do prog <- getProgName
                          let progline = prog ++ " [OPTIONS] <BAM> <OUT>"
                          hPutStr stderr $ usageInfo progline optDescrs
                          hPutStrLn stderr errs
                                               
data Conf = Conf { confBed :: !(Maybe FilePath)
                 , confWindow :: !Pos.Offset
                 , confWindowStepM :: !(Maybe Pos.Offset)
                 , confWindowMin :: !Int
                 , confVerbose :: Bool
                 } deriving (Show)

confWindowStep :: Conf -> Pos.Offset
confWindowStep conf = fromMaybe (confWindow conf) $ confWindowStepM conf

defaultWindow :: Pos.Offset
defaultWindow = 1000

defaultWindowMin :: Int
defaultWindowMin = 0

data Arg = ArgBed { unArgBed :: !String }
         | ArgASite { unArgASite :: !String }
         | ArgWindow { unArgWindow :: !String }
         | ArgWindowStep { unArgWindowStep :: !String }
         | ArgWindowMin { unArgWindowMin :: !String }
         | ArgVerbose
         deriving (Show, Read, Eq, Ord)

argBed :: Arg -> Maybe String
argBed (ArgBed bed) = Just bed
argBed _ = Nothing

argASite :: Arg -> Maybe String
argASite (ArgASite aSite) = Just aSite
argASite _ = Nothing

argWindow :: Arg -> Maybe String
argWindow (ArgWindow win) = Just win
argWindow _ = Nothing

argWindowStep :: Arg -> Maybe String
argWindowStep (ArgWindowStep wstep) = Just wstep
argWindowStep _ = Nothing

argWindowMin :: Arg -> Maybe String
argWindowMin (ArgWindowMin wmin) = Just wmin
argWindowMin _ = Nothing

optDescrs :: [OptDescr Arg]
optDescrs = [ Option ['b'] ["bed"]         (ReqArg ArgBed "BED")          "Bed filename"
            , Option []    ["window"]      (ReqArg ArgWindow "WINDOW")    ("Window size [" ++ show (fromIntegral defaultWindow :: Int) ++ "]")
            , Option []    ["window-step"] (ReqArg ArgWindowStep "STEP")  "Step between windows [WINDOW SIZE]"
            , Option []    ["window-min"]  (ReqArg ArgWindowMin "MIN")    ("Minimum counts for reporting [" ++ show defaultWindowMin ++ "]")
            , Option ['v'] ["verbose"]     (NoArg ArgVerbose)             "Verbose"
            ]

argsToConf :: [Arg] -> Either String Conf
argsToConf = runReaderT conf
    where conf = Conf <$> 
                 findBed <*>
                 findWindow <*>
                 findWindowStep <*>
                 findWindowMin <*>
                 ReaderT (return . elem ArgVerbose)
          findBed = ReaderT $ return . listToMaybe . mapMaybe argBed
          findWindow = ReaderT $ maybe (return defaultWindow) parseInt . listToMaybe . mapMaybe argWindow
          findWindowStep = ReaderT $ maybe (return Nothing) (liftM Just . parseInt)  . listToMaybe . mapMaybe argWindowStep
          findWindowMin = ReaderT $ maybe (return defaultWindowMin) parseInt . listToMaybe . mapMaybe argWindowMin
          parseInt :: (Integral a) => String -> Either String a
          parseInt = AP.parseOnly AP.decimal . BS.pack
                                               