from Bio import SeqIO

data=SeqIO.parse("milti_sequence.gb","genbank")

for s in data:
    print(s.id)
    print(s.annotations["date"])
    print(s.annotations["organism"])
    for i in s.features:
     print(i.type,"----->",i.location)