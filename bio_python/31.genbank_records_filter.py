from Bio import SeqIO

data=SeqIO.parse("insulin.gb","genbank")

rec=[]
for s in data:
    if(len(s.seq)>=1000):
        print(s)
        rec.append(s)

SeqIO.write(rec,"new_insulin.gb","genbank")