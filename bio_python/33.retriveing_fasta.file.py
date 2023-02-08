from Bio import SeqIO

data=SeqIO.parse("ls_orchid.fasta","fasta")

for s in data:

    if(len(s)<750):
        print("length less than 750: ",s.seq)
        print(">",s.id)
        print(len(s))
       






