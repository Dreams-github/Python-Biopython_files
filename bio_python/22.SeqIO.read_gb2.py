from Bio import SeqIO

data=SeqIO.read("sequence.gb","genbank")
#print(data)

print(data.id)
print(data.seq)
print(data.features)
print(data.name)
print(data.description)
