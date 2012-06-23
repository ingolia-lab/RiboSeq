module Bio.RiboSeq.Counting
       where

import Control.Applicative

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

indexEnumCount :: (Show e, Enum e) => EnumCount e -> e -> Int
indexEnumCount (EnumCount lb _ub ct) e | idx >= 0 && idx < U.length ct = ct U.! idx
                                       | otherwise = error . unwords $
                                                     [ "indexEnumCount: out of bounds: "
                                                     , show idx, show e
                                                     ]
  where idx = fromEnum e - fromEnum lb
        
enumCountTotal :: EnumCount e -> Int
enumCountTotal (EnumCount _lb _ub ct) = U.sum ct
