from Bio import Entrez
Entrez.email="swapna.2010aedhu@gmail.com"
f=open("id.txt","r")
r=open("new_id.fasta","w")
id=[]
for i in f:
    id.append(i.strip())

print(id)

for k in id:
    record=Entrez.efetch(db="nucleotide",id=k,rettype="fasta")
    print(record.read(),file=r)