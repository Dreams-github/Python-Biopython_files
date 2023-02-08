from Bio import SeqIO

data=SeqIO.parse("test.gb","genbank")

for s in data:

    for i in s.features:
        if(i.type=="CDS"):
            f=open(s.id+".fasta","w")
            s1=i.location._start.position
            s2=i.location._end.position
            print(s1,s2)
            print(s.id)
            print(s.seq[s1:s2])
            header=">"+s.id+s.description
            out=header+"\n"+s.seq[s1:s2]
            print(out,file=f)
            print(s1, s2,file=f)