from Bio import SeqIO

data=SeqIO.parse("ls_orchid.fasta","fasta")

for s in data:
    if(len(s)<=750):
        print("seq length less than 750: ",s.seq)
        id=s.id.split("|")
        print(id)
        f=open(id[-1]+".fasta","w")
        out=">"+s.id+"\n"+s.seq
        print(out,file=f)


