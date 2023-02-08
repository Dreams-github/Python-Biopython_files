f=open("sequence(3).fasta","r")

pattern=str(input("enter the pattern"))
cnt=0
for s in f:
    if(">" in s):
        print(s)
    else:
        #print(s)

        cnt=cnt+s.count(pattern)

print("count of",pattern,cnt)