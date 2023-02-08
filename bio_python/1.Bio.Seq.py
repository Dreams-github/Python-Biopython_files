from Bio import Seq
#Bio --> package   Seq--> program seq.py
from Bio.Seq import Seq
#Bio--> package .Seq--> program Seq -->class present in program

s1=Seq("ATGGTTCGCGCGATTACAGCAGCTTAAGCGTCAGAACTGAATCCAACTAGCACTGACACCT")

print(s1)
print(type(s1))


print(s1.count("A"))

print(s1.count("T"))

print(s)