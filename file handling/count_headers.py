f=open("ls_orchid.fasta","r")
c=0
for s in f:
    if(">" in s):
        print(s)
        c=c+1

print("count of headers",c)



