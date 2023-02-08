from Bio.Seq import Seq
from Bio import SeqIO

data=SeqIO.read("sequence (7).fasta","fasta")
#print(data)

print("1st frame: ",data.seq.translate())
print("2nd frame: ",data.seq[1:].translate())
print("3rd frame: ",data.seq[2:].translate())

rev=data.seq.reverse_complement()
print("4th frame: ",rev.translate())
print("5th frame: ",rev[1:].translate())
print("6th frame: ",rev[2:].translate())