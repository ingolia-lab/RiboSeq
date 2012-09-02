{-# LANGUAGE OverloadedStrings #-}
module Bio.RiboSeq.Translation
       where 
       
import Data.ByteString.Char8 as BS       

ntToAaAt :: BS.ByteString -> Int -> Maybe Amino
ntToAaAt sequ i | i < 0     = Nothing
                | otherwise = lookup (BS.take 3 $ BS.drop i sequ) trans_tbl

--translateSequ :: BS.ByteString -> [Maybe Amino]
--translateSequ ntsequ = 

data Amino = Ala | Arg | Asn | Asp | Cys | Gln | Glu | Gly
           | His | Ile | Leu | Lys | Met | Phe | Pro | Ser
           | Thr | Tyr | Trp | Val 
           | STP
           deriving (Show,Eq)
              
oneLetterM :: Maybe Amino -> Char
oneLetterM = maybe 'N' oneLetter

oneLetter :: Amino -> Char
oneLetter aa = case aa of
  Ala -> 'A'
  Arg -> 'R'
  Asn -> 'N'
  Asp -> 'D'
  Cys -> 'C'
  Gln -> 'Q'
  Glu -> 'E'
  Gly -> 'G'
  His -> 'H'
  Ile -> 'I'
  Leu -> 'L'
  Lys -> 'K'
  Met -> 'M'
  Phe -> 'F'
  Pro -> 'P'
  Ser -> 'S'
  Thr -> 'T'
  Tyr -> 'Y'
  Trp -> 'W'
  Val -> 'V'
  STP -> '*'

trans_tbl :: [(BS.ByteString,Amino)]
trans_tbl = [("GCT",Ala),("GCC",Ala),("GCA",Ala),("GCG",Ala),
             ("CGT",Arg),("CGC",Arg),("CGA",Arg),("CGG",Arg),("AGA",Arg),("AGG",Arg),
             ("AAT",Asn),("AAC",Asn),
             ("GAT",Asp),("GAC",Asp),
             ("TGT",Cys),("TGC",Cys),
             ("CAA",Gln),("CAG",Gln),
             ("GAA",Glu),("GAG",Glu),
             ("GGT",Gly),("GGC",Gly),("GGA",Gly),("GGG",Gly),
             ("CAT",His),("CAC",His),
             ("ATT",Ile),("ATC",Ile),("ATA",Ile),
             ("TTA",Leu),("TTG",Leu),("CTT",Leu),("CTC",Leu),("CTA",Leu),("CTG",Leu),
             ("AAA",Lys),("AAG",Lys),
             ("ATG",Met),
             ("TTT",Phe),("TTC",Phe),
             ("CCT",Pro),("CCC",Pro),("CCA",Pro),("CCG",Pro),
             ("TCT",Ser),("TCC",Ser),("TCA",Ser),("TCG",Ser),("AGT",Ser),("AGC",Ser),
             ("ACT",Thr),("ACC",Thr),("ACA",Thr),("ACG",Thr),
             ("TAT",Tyr),("TAC",Tyr),
             ("TGG",Trp),
             ("GTT",Val),("GTC",Val),("GTA",Val),("GTG",Val),
             ("TAG",STP),("TGA",STP),("TAA",STP)
            ]