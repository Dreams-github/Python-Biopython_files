from Bio import Entrez

Entrez.email="swapna.2010aedhu@gmail.com"

f=open("id.txt","r")


id=[]
for i in f:
    id.append(i.strip())

print(id)

for k in id:
    record=Entrez.efetch(db="nucleotide",id=k,rettype="gb")
    rec = open(k+".genbank", "w")
    print(record.read(),file=rec)
