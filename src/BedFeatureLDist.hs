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
import Statistics.Quantile

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
                                in do lp <- countTrxReads hidx trx
                                      hPutStrLn hout $ compareLine cf name lp
          in IterIO.fileDriver cfIter bed
        doCompareBam _cf _cfout = undefined
--        doCompareBam cf cfout = do ct <- countAllReads bamfile
--                                   writeFile cfout $ unlines [ compareLine cf bamfile ct ]

compareLine :: U.Vector Int -> String -> LenProf -> String
compareLine cf trx lp = intercalate "\t" [ trx
                                         , showFFloat (Just 4) score ""
                                         , show $ lpTotal lp
                                         , showFFloat (Just 4) scorel2 ""
                                         , showFFloat (Just 2) v95 ""
                                         , showFFloat (Just 1) winsTotal ""
                                         , showFFloat (Just 4) winsScore ""
                                         , showFFloat (Just 4) winsRatio ""
                                         ]
  where cfdist = U.map fromIntegral cf
        lpdist = U.map fromIntegral . lpFlatten $ lp
        score = distribL1 cfdist lpdist
        scorel2 = distribL2 cfdist lpdist
        v95 = continuousBy medianUnbiased 19  20 xprof
        v98 = continuousBy medianUnbiased 49  50 xprof
        v99 = continuousBy medianUnbiased 99 100 xprof
        xprof = case lpProfile lp of
          lprof | V.null lprof -> V.singleton 0.0           
                | otherwise -> V.map fromIntegral lprof
        wdist = lpWinsorDist lp (max v95 1.0)
        winsTotal = U.sum wdist
        winsScore = distribL1 cfdist wdist
        winsRatio = winsTotal / (fromIntegral $ lpTotal lp)

countAllReads :: FilePath -> IO (U.Vector Int)
countAllReads bamfile = do ct <- newCount
                           IterLL.run $ IterLL.joinIM $ BamIter.enumBam bamfile $ countAll ct
                           U.freeze ct
  where countAll ct = IterLL.mapM_ $ \bam -> case Bam.refSpLoc bam of
          Just loc -> countOne ct (fromIntegral $ Loc.length loc)
          Nothing -> return ()

countBedReads :: FilePath -> FilePath -> IO (U.Vector Int)
countBedReads bedfile bamfile = BamIdx.withIndex bamfile $ \hidx -> 
  let addTrx ct0 trx = do lp <- countTrxReads hidx trx
                          return $! U.zipWith (+) ct0 (lpFlatten lp)
      trxIter = Bed.bedTranscriptEnum $ IterLL.foldM addTrx (U.replicate (initlen+1) 0) 
  in IterIO.fileDriver trxIter bedfile

countTrxReads :: BamIdx.IdxHandle -> Transcript -> IO LenProf
countTrxReads hidx trx = maybe noTarget withTarget $ Bam.lookupTarget header $ unSeqLabel trxref
  where header = BamIdx.idxHeader hidx
        (OnSeq trxref trxsploc) = location trx
        trxbnds = (fromIntegral *** fromIntegral) $ Loc.bounds trxsploc
        noTarget = do when verbose $ hPutStrLn stderr . unwords $ [ "Skipping", show . unSeqLabel . geneId $ trx
                                                                  , "on", show . unSeqLabel $ trxref ]
                      return $ LenProf V.empty
        withTarget tid = do when verbose $ hPutStrLn stderr . unwords $ [ "Running", show . unSeqLabel . geneId $ trx ]
                            lprof <- newLenProf (fromIntegral . Loc.length $ trxsploc)
                            let countReads = IterLL.mapM_ $ countRead' lprof trxsploc
                                targetIter = IterLL.joinIM $ BamIter.enumIndexRegion hidx tid trxbnds countReads
                            IterLL.run targetIter
                            lpFreeze lprof
                     
countRead :: UM.IOVector Int -> SpLoc.SpliceLoc -> Bam.Bam1 -> IO ()
countRead ct trxsploc b = case Bam.refSpLoc b >>= flip SpLoc.locInto trxsploc of
  Just into | Loc.strand into == Plus -> countOne ct (fromIntegral $ Loc.length into)
  _ -> return ()

countRead' :: (MonadIO m) => LenProfIO -> SpLoc.SpliceLoc -> Bam.Bam1 -> m ()
countRead' lprof trxsploc b = case Bam.refSpLoc b >>= flip SpLoc.locInto trxsploc of
  Just into -> case Loc.startPos into of
    (Pos.Pos off Plus) -> countOne' lprof off (fromIntegral $ Loc.length into)
    _ -> return ()
  _ -> return ()

countOne' :: (MonadIO m) => LenProfIO -> Pos.Offset -> Int -> m ()
countOne' (LenProf v) off len = liftIO $ do n0 <- UM.read vct len'
                                            UM.write vct len' $! succ n0
  where len' = min len (UM.length vct - 1)
        vct = v V.! off'
        off' = max 0 $ min (fromIntegral off) (V.length v - 1)

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
          
distribL1 :: U.Vector Double -> U.Vector Double -> Double
distribL1 d1 d2 = (* 0.5) $ U.sum $ U.zipWith absDiff growD1 growD2
  where maxlen = max (U.length d1) (U.length d2)
        growD1 = (toDistrib d1) U.++ (U.replicate (maxlen - U.length d1) 0.0)
        growD2 = (toDistrib d2) U.++ (U.replicate (maxlen - U.length d2) 0.0)
        toDistrib d = U.map (/ total) d
          where total = U.sum d
        absDiff x1 x2 = abs $ x1 - x2

distribL2 :: U.Vector Double -> U.Vector Double -> Double
distribL2 d1 d2 = sqrt . U.sum . U.map (^ 2) $ U.zipWith (-) growD1 growD2
  where maxlen = max (U.length d1) (U.length d2)
        growD1 = (toDistrib d1) U.++ (U.replicate (maxlen - U.length d1) 0.0)
        growD2 = (toDistrib d2) U.++ (U.replicate (maxlen - U.length d2) 0.0)
        toDistrib d = U.map (/ total) d
          where total = U.sum d

newtype LenProfBase v e = LenProf { unLenProf :: V.Vector (v e) }
type LenProf = LenProfBase U.Vector Int
type LenProfIO = LenProfBase UM.IOVector Int

initlen :: Int
initlen = 51

zerodist :: (U.Unbox a, Num a) => U.Vector a
zerodist = U.replicate (initlen + 1) 0

newLenProf :: (MonadIO m) => Int -> m LenProfIO
newLenProf l = liftIO . liftM LenProf . V.replicateM l . UM.replicate (initlen + 1) $ 0

lpFreeze :: (MonadIO m) => LenProfIO -> m LenProf
lpFreeze = liftM LenProf . liftIO . V.mapM U.freeze . unLenProf

lpTotal :: LenProf -> Int
lpTotal = V.sum . V.map U.sum . unLenProf

lpFlatten :: LenProf -> U.Vector Int
lpFlatten (LenProf v) = V.foldl' (U.zipWith (+)) zerodist v

lpProfile :: LenProf -> V.Vector Int
lpProfile = V.map U.sum . unLenProf

lpWinsorDist :: LenProf -> Double -> U.Vector Double
lpWinsorDist (LenProf v) w = V.foldl' (U.zipWith (+)) zerodist 
                             . V.map winsor $ v
  where winsor zdist = let xdist = U.map fromIntegral zdist
                       in case U.sum xdist of
                         xttl | xttl > w -> U.map (* (w / xttl)) xdist
                              | otherwise -> xdist

data Conf = Conf { confDistOut :: !(Maybe FilePath)
                 , confCompare :: !(Maybe (FilePath, FilePath))
                 , confBedFile :: !(Maybe FilePath)
                 , confJackpot :: !(Maybe Double)
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
         | ArgJackpot { unArgJackpot :: !String }
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

argJackpot :: Arg -> Maybe String
argJackpot (ArgJackpot jp) = Just jp
argJackpot _ = Nothing

optDescrs :: [OptDescr Arg]
optDescrs = [ Option ['o'] ["output"]     (ReqArg ArgDistOut "OUTFILE")       "Length distribution output filename"
            , Option ['c'] ["compare"]    (ReqArg ArgCompareDist "DISTFILE")  "Length distribution comparison filename"
            , Option ['b'] ["bed"]        (ReqArg ArgBedFile "BEDFILE")       "Transcript bed filename"
            , Option ['s'] ["scores"]     (ReqArg ArgCompareScore "SCOREFILE") "Per-transcript comparison output filename"
            , Option ['j'] ["jackpot"]    (ReqArg ArgJackpot "Q")             "Read count quantile for jackpot suppression"
            ]

argsToConf :: [Arg] -> Either String Conf
argsToConf = runReaderT conf
    where conf = Conf <$> 
                 findDistOut <*>
                 findCompare <*>
                 findBed <*>
                 findJackpot
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
          findJackpot = ReaderT $ maybe (return Nothing) (liftM Just . parseJackpot) . listToMaybe . mapMaybe argJackpot
          parseJackpot = AP.parseOnly AP.double . BS.pack

verbose :: Bool
verbose = False