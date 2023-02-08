f=open("sequence .fasta","r")

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
        cT=cT+s.count("T")
        cG=cG+s.count("G")
        cC=cC+s.count("C")

print("count A",cA)
print("count T",cT)
print("count G",cG)
print("count C",cC)