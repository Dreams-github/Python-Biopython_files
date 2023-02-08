s="ATCTGTGCTTCGTCTGCTGCTAGTCGATTCTTCATGTGTAG"

cA=0
cT=0
cG=0
cC=0
for i in s:
    if(i=="A"):
        cA=cA+1
    elif(i=="T"):
        cT=cT+1
    elif(i=="G"):
        cG=cG+1
    elif(i=="C"):
        cC=cC+1

print("count A",cA)
print("count T",cT)
print("count G",cG)
print("count C",cC)

#count A T G C using for loop GC % AT %