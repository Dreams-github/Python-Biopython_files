s1="ATGGGGGCCAAAATTTGGGCCCATTTGGGCCCATATAAAATTGGGCCCCATTTGGGCC"
s2="CTTGAGGACAAGATGTGAGCCAATCTGTGCACACATAGAAG"


cM=0
cN=0
for i in range(len(s2)):
    if(s2[i]==s1[i]):
        print("matching seq",i)
        cM=cM+1
    elif(s2[i]!=s1[i]):
        print("not matching seq",i)
        cN=cN+1
print("count of matching",cM)
print("count of not matching",cN)
