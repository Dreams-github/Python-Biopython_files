s="ATGCTGATGGTAGCGTAAAGTGACCGTAGTAVGCCGTAACTGATT"

bp=["A","T","G","C"]

for b1 in bp:
    for b2 in bp:
        for b3 in bp:
            b=b1+b2+b3

            #print(b)
            c=s.count(b)
            if(c>0):
             print("count of bp",b,c)