module Bio.RiboSeq.StartSite
       where

import System.Directory
import System.Exit
import System.FilePath

import qualified Bio.SamTools.BamIndex as BamIndex
import qualified Data.Conduit as C
import qualified Data.Conduit.List as C

import Bio.SeqLoc.Transcript

import Bio.RiboSeq.BamFile
import Bio.RiboSeq.CodonAssignment
import Bio.RiboSeq.SVMLight
import Bio.RiboSeq.StartSVM


data LtmSample  = LtmSample  { ltmASites :: !FilePath, ltmBam :: !FilePath, chxBam :: !FilePath } deriving (Read, Show)
data LtmModel = LtmModel { ltmSamples :: ![LtmSample] } deriving (Read, Show)

data StartModel = StartModel { harrModel :: !HarrModel
                             , ltmModel :: !LtmModel
                             } deriving (Read, Show)
                                        
-- Identiﬁcation of TIS Positions. A peak is deﬁned at the nucleotide
-- level on a transcript. A peak position satisﬁes the following conditions: (i) The transcript must have both LTM and CHX reads.
-- (ii) The position must have at least 10 reads from the LTM data.
-- (iii) The position must be a local maximum within seven nucleotides (4). The position must have RLTM − RCHX of at least 0.05,
-- where Rk = (Xk/Nk) × 10 (k = LTM, CHX), Xk is the number of
-- reads on that position in data k, and Nk
-- is the total number of reads
-- on that transcript in data k. Generally, a peak position is also
-- designated a TIS. However, if a peak was not detected on the ﬁrst
-- position of any AUG or near-cognate start codon but was present
-- at the ﬁrst position of a codon immediately preceding or succeeding one of these codons, the position was designated a TIS.