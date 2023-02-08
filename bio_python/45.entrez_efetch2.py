from Bio import Entrez
from Bio import SeqIO

Entrez.email="swapna.2010aedhu@gmail.com"

data=Entrez.efetch(db="nucleotide",id="202471",rettype="gb")

rec=SeqIO.read(data,"genbank")


print(rec)
SeqIO.write(rec,"ret.gb","genbank")


