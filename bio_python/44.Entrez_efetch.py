from Bio import Entrez
from Bio import SeqIO
Entrez.email="swapna.2010aedhu@gmail.com"

data=Entrez.efetch(db="protein",id="NP_001191615.1",rettype="fasta")
rec=SeqIO.read(data,"fasta")

print(rec)

SeqIO.write(rec,"ret.fasta","fasta")