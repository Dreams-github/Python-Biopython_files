from Bio.Seq import Seq
from Bio import SeqIO
from Bio import SeqUtils

ss=SeqIO.read("sequence (7).fasta","fasta")
print(ss)

print("count of A: ",ss.seq.count("A"))
print("count of T: ",ss.seq.count("T"))
print("count of G: ",ss.seq.count("G"))
print("count of C: ",ss.seq.count("C"))
print("GC %: ",SeqUtils.GC(ss.seq))
AT=100-SeqUtils.GC(ss.seq)
print("AT %: ",AT)


