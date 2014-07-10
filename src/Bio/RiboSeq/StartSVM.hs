module Bio.RiboSeq.StartSVM
       where

import Control.Monad
import qualified Data.ByteString.Char8 as BS
import Data.Maybe
import System.Directory
import System.Exit
import System.FilePath
import System.IO
import System.Posix.Temp

import Control.Monad.CatchIO

import qualified Data.Conduit as C
import qualified Data.Conduit.List as C
import qualified Data.Vector.Unboxed as U

import qualified Bio.SamTools.BamIndex as BamIndex
import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.OnSeq
import Bio.SeqLoc.Transcript

import Bio.RiboSeq.BamFile
import Bio.RiboSeq.CodonAssignment
import Bio.RiboSeq.SVMLight

data HarrSample = HarrSample { harrBam, harrASites :: !FilePath } deriving (Eq, Ord, Read, Show)

data HarrModel = HarrModel { harrAaFields :: ![[Int]]
                           , harrSamples :: ![HarrSample]
                           , harrSvmModel :: !FilePath
                           , harrMinTotal :: !Int
                           , harrBinDir :: !FilePath
                           } deriving (Read, Show)

readHarrProfiles :: HarrModel -> Transcript -> IO TrxProfSet
readHarrProfiles harrModel trx = do profs <- mapM (readHarrProfile trx) . harrSamples $ harrModel
                                    return $! TrxProfSet trx (map (\(TrxVector _ prof) -> prof) profs)

readHarrProfile :: Transcript -> HarrSample -> IO TrxProfile
readHarrProfile trx (HarrSample bam asites) = do asd <- readASiteDelta asites
                                                 prof <- BamIndex.withIndex bam $ \hidx -> 
                                                   transcriptNtProfile asd hidx trx
                                                 return $! TrxVector trx prof

defaultHarrAaFields :: [[Int]]
defaultHarrAaFields = [[-2, -1], [0], [1], [2], [3, 4], [5, 6, 7], [8, 9, 10], [11, 12, 13]]

trainHarr :: HarrModel -> TrainPosns -> [Transcript] -> IO ExitCode
trainHarr model posns trxset = trainSVM learn fields posns (harrMinTotal model) trxsrc
  where learn = defaultLearn { learnBinary = harrBinDir model </> "svm_learn"
                             , modelOut = harrSvmModel model
                             , softMargin = Just 2.0
                             , positiveWeight = Just 4.0
                             , removeRetrain = Just True
                             , kernel = Just $ SVMRadial (Just 2.4)
                             }
        fields = aaSvmFields . harrAaFields $ model
        trxsrc = C.sourceList trxset C.$= C.mapM (readHarrProfiles model)

testHarr :: HarrModel -> TrainPosns -> Bool -> [Transcript] -> IO TestScore
testHarr model posns keeptemp trxset = testSVM classify fields posns (harrMinTotal model) keeptemp trxsrc
  where classify = defaultClassify { classifyBinary = harrBinDir model </> "svm_classify"
                                   , modelIn = harrSvmModel model
                                   }
        fields = aaSvmFields . harrAaFields $ model
        trxsrc = C.sourceList trxset C.$= C.mapM (readHarrProfiles model)

svmStartProfile :: HarrModel -> Transcript -> IO (U.Vector Bool)
svmStartProfile model trx = do profset <- readHarrProfiles model trx
                               scoreprof <- scoreProfile classify fields (harrMinTotal model) profset "."
                               return $! U.map (> 0.0) scoreprof
  where classify = defaultClassify { classifyBinary = harrBinDir model </> "svm_classify"
                                   , modelIn = harrSvmModel model
                                   }
        fields = aaSvmFields . harrAaFields $ model

scoreProfile :: RunSVMClassify -> SVMFields -> Int -> TrxProfSet -> FilePath -> IO (U.Vector Double)
scoreProfile classify fields mintotal profset scoredir
    = withTempFileName False (scoredir </> "query") $ \query ->
      withTempFileName False (scoredir </> "score") $ \rawscore -> do
        BS.writeFile query . BS.unlines $ profileQuerySet fields mintotal profset
        err <- runSVMClassify $ classify { queriesIn = query, decisionsOut = rawscore }
        case err of
          ExitSuccess -> readScoreVector query rawscore
          ExitFailure _ -> error "Failure running SVM classification"
  where readScoreVector :: FilePath -> FilePath -> IO (U.Vector Double)
        readScoreVector queryfile scorefile = do 
            query <- liftM BS.lines . BS.readFile $ queryfile
            score <- liftM BS.lines . BS.readFile $ scorefile
            return $ U.unfoldrN (length query) merge (query, score)
        merge ([], []) = Nothing
        merge ([], _rest) = error $ "scoreProfileNew: leftover scores"
        merge ((q0:qrest), [])
          | isComment q0 = Just (0/0, (qrest, []))
          | otherwise = error $ "scoreProfileNew: too few scores"
        merge ((q0:qrest), scs@(sc0:screst))
          | isComment q0 = Just (0/0, (qrest, scs))
          | otherwise = case reads . BS.unpack $ sc0 of
              [(x, "")] -> Just (x, (qrest, screst))
              _ -> error $ "scoreProfileNew: unparseable score " ++ show sc0
        isComment = maybe blankQuery ((== '#') . fst) . BS.uncons
        blankQuery = error $ "scoreProfileNew: blank query"
                
profileQuerySet :: SVMFields -> Int -> TrxProfSet -> [BS.ByteString]
profileQuerySet fields mintotal profset = map exampleAt posns
  where posns = [0..(fromIntegral . Loc.length . unOnSeq . location . transcript $ profset)]
        exampleAt p = fromMaybe justInfo . profilesExample fields mintotal targetQuery profset $ p
          where justInfo = '#' `BS.cons` profilesInfo fields profset p

trainSVM :: RunSVMLearn -> SVMFields -> TrainPosns -> Int -> C.Source IO TrxProfSet -> IO ExitCode
trainSVM learn fields posns mintotal trxsrc = withTempFile False exampleBase $ \(exfile, hex) -> do
  putTrainingSet hex fields posns mintotal trxsrc
  runSVMLearn $ learn { examplesIn = exfile }
    where exampleBase = (takeBaseName . modelOut $ learn) ++ "-training"
    
testSVM :: RunSVMClassify -> SVMFields -> TrainPosns -> Int -> Bool -> C.Source IO TrxProfSet -> IO TestScore
testSVM classify fields posns mintotal keeptemp trxsrc
    = withTempFile keeptemp exampleBase $ \(exfile, hex) ->
      withTempFileName keeptemp scoreBase $ \score ->
          do putTrainingSet hex fields posns mintotal trxsrc
             hClose hex
             ec <- runSVMClassify $ classify { queriesIn = exfile, decisionsOut = score }
             case ec of
               ExitSuccess   -> partitionTest exfile score
               ExitFailure _ -> do hPutStrLn stderr $ "SVM classification failed"
                                   return (TestData U.empty U.empty U.empty U.empty)
    where exampleBase = "svm-test-examples"
          scoreBase = "svm-test-scores"

putTrainingSet :: Handle -> SVMFields -> TrainPosns -> Int -> C.Source IO TrxProfSet -> IO ()
putTrainingSet h fields posns mintotal trxsrc = trxsrc C.$$ C.mapM_ putTraining
  where putTraining = BS.hPutStr h . BS.unlines . profileCdsExampleSet fields posns mintotal
    
profileCdsExampleSet :: SVMFields -> TrainPosns -> Int -> TrxProfSet -> [BS.ByteString]
profileCdsExampleSet fields posns mintotal profset 
  = maybe [] (cdsExampleSet fields posns mintotal profset) $
    trxProfSetStart profset

cdsExampleSet  :: SVMFields -> TrainPosns -> Int -> TrxProfSet -> Int -> [BS.ByteString]
cdsExampleSet fields (TrainPosns starts nonstarts) mintotal profs cdsstart 
  = mapMaybe (examples targetYes) starts ++ mapMaybe (examples targetNo) nonstarts
  where examples target ntoff = profilesExample fields mintotal target profs (cdsstart + ntoff)

data TrainPosns = TrainPosns { posnStarts, posnNonStarts :: ![Int] } deriving (Read, Show)

defaultTrainPosns :: TrainPosns
defaultTrainPosns = TrainPosns { posnStarts = [ 0 ]
                               , posnNonStarts = [ -6, -3, 3, 9, 18, 30, 60, 90, 120, 150 ]
                               }

profilesExample :: SVMFields -> Int -> BS.ByteString -> TrxProfSet -> Int -> Maybe BS.ByteString
profilesExample fields mintotal target profset ntpos = do 
  fvals <- profilesFields fields profset ntpos
  case sum fvals of
    ttl | ttl >= mintotal -> return $ exampleLine target (zip [1..] $ normalize fvals) (Just $ profilesInfo fields profset ntpos)
        | otherwise -> Just $! '#' `BS.cons` profilesInfo fields profset ntpos

profilesInfo :: SVMFields -> TrxProfSet -> Int -> BS.ByteString
profilesInfo fields profset ntpos = BS.unwords [ trxProfSetName profset, BS.pack . show $ ntpos, readCountBS, groupBS, posnsBS ]
    where readCountBS = maybe (BS.pack "n/a") (BS.pack . show . sum) $ profilesFields fields profset ntpos
          fieldsBS = BS.pack . show $ profilesFields fields profset ntpos
          groupBS = BS.pack . show $ map (\fgr -> profilesFields fgr profset ntpos) $ groupFields fields
          posnsBS = BS.pack . show $ profilesFields (singleFields fields) profset ntpos

profilesFields :: SVMFields -> TrxProfSet -> Int -> Maybe [Int]
profilesFields fields profs ntpos = liftM concat . mapM (vectorFields fields ntpos) $ profiles profs

data SVMFields = SVMFields { ntoffsets :: [[Int]] } deriving (Show)

minField, maxField :: SVMFields -> Int
minField = minimum . minimum . ntoffsets
maxField = maximum . maximum . ntoffsets

singleFields :: SVMFields -> SVMFields
singleFields fields = SVMFields . map (: []) $ [(minField fields)..(maxField fields)]

groupFields :: SVMFields -> [SVMFields]
groupFields fields = map (SVMFields . map (: [])) . ntoffsets $ fields

profilesTotal :: SVMFields -> TrxProfSet -> Int -> Maybe Int
profilesTotal fields profs = liftM sum . profilesFields fields profs

aaSvmFields :: [[Int]] -> SVMFields
aaSvmFields aaoffsets = SVMFields { ntoffsets = map tont aaoffsets }
    where tont = concatMap (\aaoff -> let ntctr = 3 * aaoff in [ ntctr - 1, ntctr, ntctr + 1 ])

lowestoffset :: SVMFields -> Int
lowestoffset = minimum . map minimum . ntoffsets

highestoffset :: SVMFields -> Int
highestoffset = maximum . map maximum . ntoffsets

vectorFields :: (U.Unbox a, Num a) => SVMFields -> Int -> U.Vector a -> Maybe [a]
vectorFields fields ntpos vec | inflanking = Just fieldvalues
                              | otherwise = Nothing
    where inflanking = (lowestoffset fields + ntpos >= 0)
                       && (highestoffset fields + ntpos < U.length vec)
          fieldvalues = map fieldvalue $ ntoffsets fields
          fieldvalue = sum . map (countat . (+ ntpos))
          countat ntoff | ntoff < 0 || ntoff >= U.length vec = error $ "vectorFields: " ++ show (ntoff, bounds)
                        | otherwise = vec U.! ntoff
          bounds = ( ( 0 :: Int, U.length vec - 1), (lowestoffset fields, highestoffset fields), ntpos)

normalize :: (Integral a) => [a] -> [Double]
normalize zs | total == 0 = replicate (length zs) 0.0
             | otherwise = map norm zs
    where total = sum . map (^ (2 :: Int)) $ zs
          norm = (/ (sqrt . fromIntegral $ total)) . fromIntegral

withTempFile :: Bool -> FilePath -> ((FilePath, Handle) -> IO a) -> IO a
withTempFile keepTemp tmpbase = bracket (mkstemp tmpnam) closeAndRemove 
    where tmpnam = tmpbase ++ "XXXXXX"
          closeAndRemove (tmpfile, htmp) = do hIsOpen htmp >>= flip when (hClose htmp)
                                              unless keepTemp $ doesFileExist tmpfile >>= flip when (removeFile tmpfile)

withTempFileName :: Bool -> FilePath -> (FilePath -> IO a) -> IO a
withTempFileName keepTemp tmpbase = withTempFile keepTemp tmpbase . closeHandle
    where closeHandle act (tmpfile, htmp) = hClose htmp >> act tmpfile

--        Profile.writeVectorProfileWith showf scorefile $ Profile.fromInfo info scorevec
--        readTotalProfile fields samples kgid >>= Profile.writeVectorProfileWith show totalfile
--    where scorefile = scoredir </> trxname ++ "_score.txt"
--          totalfile = scoredir </> trxname ++ "_reads.txt"
--          trxname = BS.unpack . trxProfileName $ prof0

-- readTotalProfile :: SVMFields -> [FilePath] -> BS.ByteString -> IO (Profile.Profile (U.Vector Int))
-- readTotalProfile fields samples kgid = liftM totalProfile $ Profile.readGeneProfiles samples kgid
--     where totalProfile profs@(prof0:_) = Profile.fromInfo (Profile.toInfo prof0) . U.fromListN (len) . map totalAt $ [0..(len-1)]
--               where len = U.length . Profile.profile $ prof0
--                     totalAt = fromMaybe 0 . profilesTotal fields profs

-- writeQueryProfile :: SVMFields -> [FilePath] -> BS.ByteString -> FilePath -> IO Profile.ProfileInfo
-- writeQueryProfile fields samples kgid exampleout = do profs <- Profile.readGeneProfiles samples kgid
--                                                       let posns = [0..(proflen profs - 1)]
--                                                       withFile exampleout WriteMode $ \hout ->
--                                                           mapM_ (hPutWordsLn hout . exampleAt profs) posns
--                                                       return . Profile.toInfo . head $ profs
--     where exampleAt profs posn = fromMaybe justInfo . profilesExample fields targetQuery profs $ posn
--               where justInfo = [ BS.singleton '#', profilesInfo fields profs posn ]
--           proflen (prof0:_) = U.length . Profile.profile $ prof0
--           proflen [] = error "writeQueryProfile: Empty profile list"

-- hPutWordsLn :: Handle -> [BS.ByteString] -> IO ()
-- hPutWordsLn hout = (>> hPutChar hout '\n') . sequence_ . intersperse (hPutChar hout ' ') . map (BS.hPutStr hout)

-- trainSVM :: RunSVMLearn -> SVMFields -> [FilePath] -> [BS.ByteString] -> IO ()
-- trainSVM learn fields samples  kgids
--     = withTempFileName exampleBase $ \examples -> do
--         writeTrainingSet fields kgids samples examples
--         runSVMLearn $ learn { examplesIn = examples }
--         return ()
--     where exampleBase = (takeBaseName . modelOut $ learn) ++ "-training"

-- writeTrainingSet :: SVMFields -> [BS.ByteString] -> [FilePath] -> FilePath -> IO ()
-- writeTrainingSet fields = writeCdsExampleSet fields trainingNtPosns

-- writeCdsExampleSet :: SVMFields -> ([Int], [Int]) -> [BS.ByteString] -> [FilePath] -> FilePath -> IO ()
-- writeCdsExampleSet fields posns kgids samples outname = withFile outname WriteMode $ \hout ->
--                                                         forM_ kgids $ writeExamplesKgid hout
--     where writeExamplesKgid hout kgid = do profs <- Profile.readGeneProfiles samples kgid
--                                            maybe skip process . profileCdsExampleSet fields posns $ profs
--               where skip = BS.hPutStrLn stderr $ kgid `BS.append` (BS.pack "\tskipped")
--                     process ls = mapM_ (hPutWordsLn hout) ls >> BS.hPutStrLn stderr kgid

