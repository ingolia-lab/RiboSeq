module Bio.RiboSeq.Counting
       where

import Control.Applicative
import Control.Monad

import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM

data EnumCountBase v e = EnumCount { ecLB :: !e, ecUB :: !e, ecCount :: v }
type EnumCountIO e = EnumCountBase (UM.IOVector Int) e
type EnumCount e = EnumCountBase (U.Vector Int) e

newEnumCountIO :: (Enum e) => (e, e) -> IO (EnumCountIO e)
newEnumCountIO (lb, ub) = EnumCount lb ub <$> UM.replicate (1 + fromEnum ub - fromEnum lb) 0

newEnumCount :: (Enum e) => (e, e) -> EnumCount e
newEnumCount (lb, ub) = EnumCount lb ub (U.replicate (1 + fromEnum ub - fromEnum lb) 0)

countEnum :: (Show e, Enum e) => EnumCountIO e -> e -> IO ()
countEnum (EnumCount lb _ub ct) e | idx >= 0 && idx < UM.length ct = incr ct idx
                                  | otherwise = ioError . userError . unwords $
                                                [ "countEnum: out of bounds: "
                                                , show idx, show e
                                                ]
  where idx = fromEnum e - fromEnum lb
        incr v i = UM.read v i >>= UM.write v i . (succ $!)

freezeEnumCount :: EnumCountIO e -> IO (EnumCount e)
freezeEnumCount (EnumCount lb ub ctio) = EnumCount lb ub <$> U.freeze ctio

readEnumCount :: (Show e, Enum e) => EnumCountIO e -> e -> IO Int
readEnumCount (EnumCount lb _ub ct) e | idx >= 0 && idx < UM.length ct = UM.read ct idx
                                      | otherwise = ioError . userError . unwords $
                                                    [ "readEnumCount: out of bounds: "
                                                    , show idx, show e
                                                    ]
  where idx = fromEnum e - fromEnum lb

addEnumCount :: (Ord e, Enum e) => EnumCount e -> EnumCount e -> EnumCount e
addEnumCount (EnumCount lb0 ub0 ct0) (EnumCount lb1 ub1 ct1)
  = EnumCount lb ub $ U.generate (1 + fromEnum ub - fromEnum lb) ctsum
  where lb = min lb0 lb1
        ub = max ub0 ub1
        ctsum i = let i0 = i + fromEnum lb - fromEnum lb0
                      i1 = i + fromEnum lb - fromEnum lb1
                      n0 = if (i0 >= 0 && i0 < U.length ct0) then ct0 U.! i0 else 0
                      n1 = if (i1 >= 0 && i1 < U.length ct1) then ct1 U.! i1 else 0
                  in n0 + n1
        
indexEnumCount :: (Show e, Enum e) => EnumCount e -> e -> Int
indexEnumCount (EnumCount lb _ub ct) e | idx >= 0 && idx < U.length ct = ct U.! idx
                                       | otherwise = error . unwords $
                                                     [ "indexEnumCount: out of bounds: "
                                                     , show idx, show e
                                                     ]
  where idx = fromEnum e - fromEnum lb
        
enumCountTotal :: EnumCount e -> Int
enumCountTotal (EnumCount _lb _ub ct) = U.sum ct

newtype ProfileEnumCountBase v e = ProfileEnumCount { pecProf :: (V.Vector (EnumCountBase v e)) }
type ProfileEnumCountIO e = ProfileEnumCountBase (UM.IOVector Int) e
type ProfileEnumCount e = ProfileEnumCountBase (U.Vector Int) e

newProfileEnumCountIO :: (Enum e) => Int -> (e, e) -> IO (ProfileEnumCountIO e)
newProfileEnumCountIO len bnds = liftM ProfileEnumCount $ V.replicateM len (newEnumCountIO bnds)

profileCountEnum :: (Show e, Enum e) => ProfileEnumCountIO e -> Int -> e -> IO ()
profileCountEnum (ProfileEnumCount prof) i | i >= 0 && i < V.length prof = countEnum (prof V.! i)
                                           | otherwise = const . ioError . userError . unwords $
                                                         [ "profileCountEnum: out of profile: "
                                                         , show i, show $ V.length prof
                                                         ]
                                                         
freezeProfileEnumCount :: ProfileEnumCountIO e -> IO (ProfileEnumCount e)
freezeProfileEnumCount (ProfileEnumCount profio) = liftM ProfileEnumCount $
                                                   V.mapM freezeEnumCount profio
                                                   
indexProfileEnumCount :: (Show e, Enum e) => ProfileEnumCount e -> Int -> e -> Int
indexProfileEnumCount (ProfileEnumCount prof) i | i >= 0 && i < V.length prof = indexEnumCount (prof V.! i)
                                                | otherwise = error . unwords $
                                                              [ "indexProfileEnumCount: out of profile: "
                                                              , show i, show $ V.length prof
                                                              ]
                                                              
profileEnumCountTotal :: ProfileEnumCount e -> Int
profileEnumCountTotal (ProfileEnumCount prof) = V.sum . V.map enumCountTotal $ prof

profileEnumCountTotalProfile :: ProfileEnumCount e -> V.Vector Int
profileEnumCountTotalProfile (ProfileEnumCount prof) = V.map enumCountTotal prof

profileEnumCountTotalCount :: (Ord e, Enum e) => ProfileEnumCount e -> EnumCount e
profileEnumCountTotalCount (ProfileEnumCount prof) = V.foldl1' addEnumCount prof