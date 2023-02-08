from Bio import Entrez
code=input("enter the code")
Entrez.email="swapna.2010aedhu@gmail.com"

data=Entrez.einfo(db=code)

info=Entrez.read(data)

print(info)