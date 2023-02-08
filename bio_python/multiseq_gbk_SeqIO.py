from Bio import SeqIO

data=SeqIO.parse("ls_orchid.gbk","genbank")

for s in data:
    print("id",s.id)
    print("seq",s.seq)
    print("features",s.features)
    print("name",s.name)
    print("annotations",s.annotations)
    print("seq translate: ",s.seq.translate())
