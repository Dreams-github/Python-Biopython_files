from Bio import Entrez
Entrez.email="swapna.2010aedhu@gmail.com"

data=Entrez.esearch(db="nucleotide",term="insulin AND Mus musculus",retmax=100,idtype="acc")

rec=Entrez.read(data)

print(rec)

print(len(rec["IdList"]))