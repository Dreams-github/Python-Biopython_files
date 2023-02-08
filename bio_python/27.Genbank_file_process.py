from Bio import SeqIO

data=SeqIO.parse("milti_sequence.gb","genbank")

for s in data:
    print(s.id)
    print(s.seq)
    print(s.annotations)
    print(s.description)
    print(s.name)
    print(s.features)