from Bio import Entrez

Entrez.email="swapna.2010aedhu@gmail.com"

data=Entrez.esearch(db="protein",term="Homo sapiens insulin AND Homo sapiens",retmax=50,idtype="acc")

rec=Entrez.read(data)

print(rec)
print(rec["IdList"])