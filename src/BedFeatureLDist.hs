module Main
  where

import Control.Applicative
import Control.Arrow
import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.Reader
import qualified Data.ByteString.Char8 as BS
import qualified Data.ByteString.Lazy.Char8 as LBS
import Data.IORef
import Data.List
import Data.Maybe
import Numeric
import System.Console.GetOpt
import System.Environment
import System.FilePath
import System.IO

import qualified Data.Attoparsec.Char8 as AP
import qualified Data.Iteratee.IO as IterIO
import qualified Data.Iteratee.ListLike as IterLL
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM

import qualified Bio.SamTools.Bam as Bam
import qualified Bio.SamTools.BamIndex as BamIdx
import qualified Bio.SamTools.Iteratee as BamIter
import qualified Bio.SeqLoc.Bed as Bed
import Bio.SeqLoc.LocRepr
import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.Position as Pos
import qualified Bio.SeqLoc.SpliceLocation as SpLoc
import Bio.SeqLoc.Strand
import Bio.SeqLoc.Transcript

main :: IO ()
main = getArgs >>= handleOpt . getOpt RequireOrder optDescrs
    where handleOpt (_,    _,       errs@(_:_)) = usage (unlines errs)
          handleOpt (_,    [],      [])         = usage "Specify a BAM file"
          handleOpt (args, [bam],   [])         = either usage (doLengthDist bam) $ argsToConf args
          handleOpt (_,    (_:_:_), [])         = usage "Specify just one BAM file"
          usage errs = do prog <- getProgName
                          let progline = prog ++ " [OPTIONS] <BAM>"
                          hPutStr stderr $ usageInfo progline optDescrs
                          hPutStrLn stderr errs

doLengthDist :: FilePath -> Conf -> IO ()
doLengthDist bamfile conf = do maybe (return ()) doDist $ confDistOut conf
                               maybe (return ()) doCompare $ confCompare conf
  where doDist distout = (maybe countAllReads countBedReads $ confBedFile conf) bamfile >>= 
                         writeFile distout . countTable
        doCompare (cfin, cfout) = do 
          cf <- liftM parseCountTable $ BS.readFile cfin
          (maybe doCompareBam doCompareBed $ confBedFile conf) cf cfout
        doCompareBed bed cf cfout = BamIdx.withIndex bamfile $ \hidx -> 
          withFile cfout WriteMode $ \hout -> 
          let cfIter = Bed.bedTranscriptEnum $ IterLL.mapM_ compareDist
              compareDist trx = let name = BS.unpack . unSeqLabel . geneId $ trx
                                in do ct <- newCount
                                      countTrxReads hidx ct trx
                                      U.freeze ct >>= hPutStrLn hout . compareLine cf name
          in IterIO.fileDriver cfIter bed
        doCompareBam cf cfout = do ct <- countAllReads bamfile
                                   writeFile cfout $ unlines [ compareLine cf bamfile ct ]

compareLine :: U.Vector Int -> String -> U.Vector Int -> String
compareLine cf trx ctvec = intercalate "\t" [ trx
                                            , showFFloat (Just 4) score ""
                                            , show $ U.sum ctvec
                                            , showFFloat (Just 4) (distribL2 cf ctvec) ""
                                            ]
  where score = distribL1 cf ctvec

countAllReads :: FilePath -> IO (U.Vector Int)
countAllReads bamfile = do ct <- newCount
                           IterLL.run $ IterLL.joinIM $ BamIter.enumBam bamfile $ countAll ct
                           U.freeze ct
  where countAll ct = IterLL.mapM_ $ \bam -> case Bam.refSpLoc bam of
          Just loc -> countOne ct (fromIntegral $ Loc.length loc)
          Nothing -> return ()

countBedReads :: FilePath -> FilePath -> IO (U.Vector Int)
countBedReads bedfile bamfile = BamIdx.withIndex bamfile $ \hidx -> do
  ct <- newCount
  let trxIter = Bed.bedTranscriptEnum . IterLL.mapM_ $ countTrxReads hidx ct 
  IterIO.fileDriver trxIter bedfile
  U.freeze ct

countTrxReads :: BamIdx.IdxHandle -> UM.IOVector Int -> Transcript -> IO ()
countTrxReads hidx ct trx = maybe noTarget withTarget $ Bam.lookupTarget header $ unSeqLabel trxref
  where header = BamIdx.idxHeader hidx
        (OnSeq trxref trxsploc) = location trx
        trxbnds = (fromIntegral *** fromIntegral) $ Loc.bounds trxsploc
        noTarget = when verbose $ hPutStrLn stderr . unwords $ [ "Skipping", show . unSeqLabel . geneId $ trx
                                                               , "on", show . unSeqLabel $ trxref ]
        withTarget tid = do when verbose $ hPutStrLn stderr . unwords $ [ "Running", show . unSeqLabel . geneId $ trx ]
                            IterLL.run targetIter
          where targetIter = IterLL.joinIM $ BamIter.enumIndexRegion hidx tid trxbnds countReads
        countReads = IterLL.mapM_ $ countRead ct trxsploc
                     
countRead :: UM.IOVector Int -> SpLoc.SpliceLoc -> Bam.Bam1 -> IO ()
countRead ct trxsploc b = case Bam.refSpLoc b >>= flip SpLoc.locInto trxsploc of
  Just into | Loc.strand into == Plus -> countOne ct (fromIntegral $ Loc.length into)
  _ -> return ()

newCount :: (MonadIO m) => m (UM.IOVector Int)
newCount = liftIO $ UM.replicate (initlen + 1) 0
  where initlen = 51

countOne :: (MonadIO m) => UM.IOVector Int -> Int -> m ()
countOne ct len = liftIO $ do l0 <- UM.read ct len'
                              UM.write ct len' $! succ l0
  where len' = min len (UM.length ct - 1)

countTable :: U.Vector Int -> String
countTable ctvec = unlines . V.toList . V.imap lenline . U.convert $ ctvec
  where lenline idx n = intercalate "\t" [ show idx
                                         , show n
                                         , showFFloat (Just 4) (fromIntegral n / total) ""
                                         , showFFloat (Just 4) (fromIntegral (cumul U.! idx) / total) ""
                                         ]
        total = fromIntegral $ U.sum ctvec
        cumul = U.prescanl' (+) 0 ctvec

parseCountTable :: BS.ByteString -> U.Vector Int
parseCountTable = U.fromList . zipWith parseCount [0..] . BS.lines
  where parseCount idx l = let bad = error $ "Bad count table line " ++ show idx ++ ": " ++ show l
                           in case BS.split '\t' l of
                             (idxbs:nbs:_) -> let idx' = either bad id $! AP.parseOnly AP.decimal idxbs
                                                  n = either bad id $! AP.parseOnly AP.decimal nbs
                                              in if idx' == idx
                                                 then n
                                                 else bad
                             _ -> bad
          
distribL1 :: U.Vector Int -> U.Vector Int -> Double
distribL1 d1 d2 = (* 0.5) $ U.sum $ U.zipWith absDiff growD1 growD2
  where maxlen = max (U.length d1) (U.length d2)
        growD1 = (toDistrib d1) U.++ (U.replicate (maxlen - U.length d1) 0.0)
        growD2 = (toDistrib d2) U.++ (U.replicate (maxlen - U.length d2) 0.0)
        toDistrib d = U.map ((/ total) . fromIntegral) d
          where total = fromIntegral $ U.sum d
        absDiff x1 x2 = abs $ x1 - x2

distribL2 :: U.Vector Int -> U.Vector Int -> Double
distribL2 d1 d2 = sqrt . U.sum . U.map (^ 2) $ U.zipWith (-) growD1 growD2
  where maxlen = max (U.length d1) (U.length d2)
        growD1 = (toDistrib d1) U.++ (U.replicate (maxlen - U.length d1) 0.0)
        growD2 = (toDistrib d2) U.++ (U.replicate (maxlen - U.length d2) 0.0)
        toDistrib d = U.map ((/ total) . fromIntegral) d
          where total = fromIntegral $ U.sum d

data Conf = Conf { confDistOut :: !(Maybe FilePath)
                 , confCompare :: !(Maybe (FilePath, FilePath))
                 , confBedFile :: !(Maybe FilePath)
                 } deriving (Show)

confCompareDistFile :: Conf -> Maybe FilePath
confCompareDistFile = liftM fst . confCompare

confCompareScoreFile :: Conf -> Maybe FilePath
confCompareScoreFile = liftM snd . confCompare

confBedFileE :: Conf -> FilePath
confBedFileE conf = fromMaybe noBedFile $ confBedFile conf
  where noBedFile = error "No bed file"

data Arg = ArgDistOut { unArgDistOut :: !String }
         | ArgCompareDist { unArgCompareDist :: !String }
         | ArgCompareScore { unArgCompareScore :: !String }
         | ArgBedFile { unArgBedFile :: !String }         
         deriving (Show, Read, Eq, Ord)

argDistOut :: Arg -> Maybe String
argDistOut (ArgDistOut out) = Just out
argDistOut _ = Nothing

argCompareDist :: Arg -> Maybe String
argCompareDist (ArgCompareDist cdist) = Just cdist
argCompareDist _ = Nothing

argCompareScore :: Arg -> Maybe String
argCompareScore (ArgCompareScore cscore) = Just cscore
argCompareScore _ = Nothing

argBedFile :: Arg -> Maybe String
argBedFile (ArgBedFile bed) = Just bed
argBedFile _ = Nothing

optDescrs :: [OptDescr Arg]
optDescrs = [ Option ['o'] ["output"]     (ReqArg ArgDistOut "OUTFILE")       "Length distribution output filename"
            , Option ['c'] ["compare"]    (ReqArg ArgCompareDist "DISTFILE")  "Length distribution comparison filename"
            , Option ['b'] ["bed"]        (ReqArg ArgBedFile "BEDFILE")       "Transcript bed filename"
            , Option ['s'] ["scores"]     (ReqArg ArgCompareScore "SCOREFILE") "Per-transcript comparison output filename"
            ]

argsToConf :: [Arg] -> Either String Conf
argsToConf = runReaderT conf
    where conf = Conf <$> 
                 findDistOut <*>
                 findCompare <*>
                 findBed
          findDistOut = ReaderT $ return . listToMaybe . mapMaybe argDistOut
          findBed = ReaderT $ return . listToMaybe . mapMaybe argBedFile
          findCompareDist = ReaderT $ return . listToMaybe . mapMaybe argCompareDist
          findCompareScore = ReaderT $ return . listToMaybe . mapMaybe argCompareScore
          findCompare = do mdist <- findCompareDist
                           mscore <- findCompareScore
                           case (mdist, mscore) of
                             (Just dist, Just score) -> return $! Just (dist, score)
                             (Nothing, Nothing) -> return Nothing
                             (Just dist, Nothing) -> lift $ Left "No comparison output filename"
                             (Nothing, Just score) -> lift $ Left "No comparison distribution given"


verbose :: Bool
verbose = False