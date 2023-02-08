s="ATGCGTGTCGTAGTAAAGTGAGCCCGTAAAGTTTGGG"
base=("A","T","G","C")
d={}

for b1 in base:
    for b2 in base:
        for b3 in base:
            b=b1+b2+b3
            c=s.count(b)
            if(c>0):
                d[b]=c
print(d)
