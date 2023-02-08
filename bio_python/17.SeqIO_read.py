from Bio import SeqIO

data=SeqIO.read("sequence (6).fasta","fasta")

#print(data)

print(data.id)

print(data.seq)

print(data.name)

print(data.description)