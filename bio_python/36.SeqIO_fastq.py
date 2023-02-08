from Bio import SeqIO

data=SeqIO.parse("sample.fastq","fastq")

for s in data:
    print(s.id)
    print(s.seq)
    print(s.letter_annotations["phred_quality"])