from Bio import Entrez

Entrez.email="swapna.2010aedhu@gmial.com"

data=Entrez.esummary(db="nucleotide",id="M57671.1,NM_001204686")

rec=Entrez.read(data)

print(rec)
f=open("summary.txt","w")
for i in rec:
    print("-----------------------------")
    for k,v in i.items():
        print(k,":",v,file=f)

