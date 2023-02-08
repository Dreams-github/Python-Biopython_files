from Bio import SeqIO
from matplotlib import pylab
from pylab import *

data=SeqIO.parse("ls_orchid.fasta","fasta")
gc=[]
for i in data:
    gc.append(len(i.seq))
#print(gc)
gc.sort()
print(gc)


pylab.hist(gc, bins=20)
pylab.title(
    "%i orchid sequences\nLengths %i to %i" % (len(gc), min(gc), max(gc))
)
pylab.xlabel("Sequence length (bp)")
pylab.ylabel("Count")
pylab.show()