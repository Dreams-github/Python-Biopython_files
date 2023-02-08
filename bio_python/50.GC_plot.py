from Bio import SeqIO
from Bio import SeqUtils
from matplotlib import pylab
from pylab import *

data=SeqIO.parse("ls_orchid.fasta","fasta")
gc=[]
for i in data:
    gc.append(SeqUtils.GC(i.seq))
#print(gc)
gc.sort()
print(gc)

pylab.plot(gc)
pylab.title("GC plot")
pylab.xlabel("Genes")
pylab.ylabel("GC %")
pylab.show()
