from Bio import SeqIO

data=SeqIO.parse("sample.fastq","fastq")

c=0
for s in data:
    #print(s.id)
    c=c+1

print("no of reads",c)