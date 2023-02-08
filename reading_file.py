f=open("C:/Users/avksa/PycharmProjects/pythonProject/file handling/sequence(3).fasta","r")

cA=0
cT=0
cG=0
cC=0
for s in f:
    if(">" in s):
        print(s)
    else:
        print(s)
        cA=cA+s.count("A")
        cT=cT + s.count("T")
        cG = cG + s.count("G")
        cC = cC + s.count("C")
print("count of A",cA)
print("count of T",cT)
print("count of G",cG)
print("count of C",cC)
l=cA+cT+cG+cC
gc=cG+cC
at=cA+cT
PGC=(gc/l)*100
PAT=(at/l)*100
print("GC %",PGC)
print("AT %",PAT)

