from Bio import SeqIO

data=SeqIO.parse("insulin.gb","genbank")

for s in data:
    for i in s.features:
        if(i.type=="CDS"):
            id="".join(i.qualifiers.get("protein_id"))
            seq="".join(i.qualifiers.get("translation"))
            f=open(id+".fasta","w")
            out=">"+id+"\n"+seq
            print(out,file=f)


#filter fasta file for length seq >300