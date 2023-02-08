from Bio.Seq import Seq
from Bio import SeqUtils

data=Seq("TCGATCGTCTAGTCGTCAGTCGTACGTCGTAGCTGACTGACGTCGTAGCTGCTGATGCTGTACGCTGATGCTGATGTAGCTGCAGCT")

print(data)

print("GC %",SeqUtils.GC(data))

print("molecular weight",SeqUtils.molecular_weight(data))

AT=100-SeqUtils.GC(data)

print("AT %",AT)