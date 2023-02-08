from Bio import SeqIO

data=SeqIO.parse("sequence.multi.fasta","fasta")

for s in data:
    print(s.id)
    print(s.seq)
    print(s.description)
    print(s.name)