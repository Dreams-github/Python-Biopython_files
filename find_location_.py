s="CATATGGTACGGTATGGACGATGGATCGATGGATGCATGGATGCTCGAATGGTACTCAATGGATCGATGGATCGTGCAATGGATCGCAATGGACTGTACATGG"
#  0123456
pattern="ATGG"
c=0
for i in range(len(s)):
    if(s.startswith(pattern,i)):
        print("start location of",pattern,i,"end location of",pattern,i+len(pattern)-1)
        c=c+1

print("count of",pattern,c)