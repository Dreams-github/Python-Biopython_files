seq="""AGTTCTAGTTCTATCTTTTTCCTGCATCCTGTCTGGAAGTTAGAAGGAAACAGACCACAGACCTGGTCCC
CAAAAGAAATGGAGGCAATAGGTTTTGAGGGGCATGAGGACGGGGTTCAGCCTCCAGGGTCCTACACACA
AATCAGTCAGTGGCCCAGAAGACCCCCCTCGGAATCGGAGCAGGGAGGATGGGGAGTGTGAGGGGTATCC
TTGATGCTTGTGTGTCCCCAACTTTCCAAATCCCCGCCCCCGCGATGGAGAAGAAACCGAGACAGAAGGT
GCAGGGCCCACTACCGCTTCCTCCAGATGAGCTCATGGG"""

print(seq)
print(len(seq))
s=seq.replace("\n","")
print(s)
print(len(s))

cA=0
cT=0
cG=0
cC=0

for j in seq:
    if(j=="A"):
        cA=cA+1
    elif(j=="T"):
        cT=cT+1
    elif(j=="G"):
        cG=cG+1
    elif(j=="C"):
        cC=cC+1
print("count of A",cA)
print("count of T",cT)
print("count of G",cG)
print("count of C",cC)

PGC=((cG+cC)/len(seq))*100
print("GC % is",PGC)

PAT=((cA+cT)/len(seq))*100
print("AT % is",PAT)

ctA=0
ctT=0
ctG=0
ctC=0
for i in s:
    if(i=="A"):
        ctA=ctA+1
    elif(i=="T"):
        ctT=ctT+1
    elif(i=="G"):
        ctG=ctG+1
    elif(i=="C"):
        ctC=ctC+1
print("count of A",ctA)
print("count of T",ctT)
print("count of G",ctG)
print("count of C",ctC)

gc=(ctG+ctC)
at=(ctA+ctT)

perGC=(gc/len(s))*100
print("GC % is",perGC)

perAT=(at/len(s))*100
print("AT % is",perAT)




