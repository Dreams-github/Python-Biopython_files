s="ATGACTCTGATGACTGATCGCATAGCTGACTGCTGCATGCATAGCTGACGTACGTCGATCGATGCTCGTAGCTC"

base=["A","T","G","C"]

for b1 in base:
    for b2 in base:
        b=b1+b2
        #print(b)
        c=s.count(b)
        print("count of",b,c)