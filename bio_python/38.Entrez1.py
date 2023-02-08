from Bio import Entrez

Entrez.email="swapna.2010aedhu@gmail.com"

data=Entrez.einfo()

info=Entrez.read(data)

print(info)