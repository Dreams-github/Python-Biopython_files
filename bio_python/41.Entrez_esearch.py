from Bio import Entrez

Entrez.email="swapna.2010aedhu@gmail.com"

data=Entrez.esearch(db="nucleotide",term="Homo sapiens insulin AND Homo sapiens",retmax=50)

rec=Entrez.read(data)

print(rec)

print(rec["IdList"])