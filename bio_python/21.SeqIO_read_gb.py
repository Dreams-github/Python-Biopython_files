from Bio import SeqIO

data=SeqIO.read("seq.gb","genbank")

#print(data)

print(data.id)

print(data.seq)

print(data.description)

print(data.features)