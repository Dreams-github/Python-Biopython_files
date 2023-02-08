from Bio.Seq import Seq
from Bio import SeqUtils
from Bio import SeqIO

data=SeqIO.parse("sequence.multi.fasta","fasta")

for s in data:
    print(s.id)
    print(s.name)
    print(s.seq)
    print(s.description)
    print(s.seq.count("A"))
    a=s.seq.count("A")
    print(s.seq.count("T"))
    t=s.seq.count("T")
    print(s.seq.count("G"))
    g=s.seq.count("G")
    print(s.seq.count("C"))
    c=s.seq.count("C")
    print("GC %:",SeqUtils.GC(s.seq))
    AT=100-SeqUtils.GC(s.seq)
    print("AT %:",AT)
    print("A%:",a/len(s.seq)*100)
    print("T%:",t/len(s.seq)*100)
    print("G%:",g/len(s.seq)*100)
    print("C%:",c/len(s.seq)*100)
