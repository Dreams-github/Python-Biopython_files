from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

data=SeqRecord(Seq("ATGGGCCTGGTCCAATTGGTCCAATGGTG"),id="NC87887",name="abcdde",description="dna seq data")

#print(data)

print(data.id)
print(data.seq)
print(data.name)
print(data.descripton)
