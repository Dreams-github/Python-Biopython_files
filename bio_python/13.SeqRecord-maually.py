from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sdata=SeqRecord(Seq("TGACGTCGATCGGCTAGGATCGAGCGCTTTACGGGACTAGGCGATGTGGCAAGCTGAGTGG"),id="G786453",name="GENE",description="Seq of some Vertebrate")
#print(sdata)

print(sdata.id)
print(sdata.name)
print(sdata.seq)
print(sdata.description)

