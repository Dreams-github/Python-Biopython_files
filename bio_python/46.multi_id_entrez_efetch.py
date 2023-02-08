from Bio import Entrez
from Bio import SeqIO

f=open("multi.fasta","w")
Entrez.email="swapna.2010aedhu@gmail.com"

data=Entrez.efetch(db="nucleotide",id="U03610.1,M57671.1",rettype="fasta")
print(data.read(),file=f)
rec=SeqIO.parse(data,"fasta")
print(rec)

for i in rec:
    print(i)


