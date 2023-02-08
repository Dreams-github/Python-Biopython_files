s="ATGGCCGGTTTAACCGGTTTGGCCAAATTGGGCCCC"

base=["A","T","G","C"]
d={}
for b1 in base:
    for b2 in base:
        b=b1+b2
        c=s.count(b)
        d[b]=c

print(d)
