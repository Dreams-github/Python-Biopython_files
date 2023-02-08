from Bio import SeqIO

data=SeqIO.read("sequence (7).fasta","fasta")

#print(data)

print(data.id)
print(data.description)