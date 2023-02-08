from Bio.Seq import Seq

s=Seq("ATGGCCAACCGGTTTGGTTGGCCAAATTGGCCCGGTTTTAAACCCGGTTTGGCCAAAACCCGGTTTAAAAA")

print(s)

print("1st frame: ",s.translate())

print("2nd frame: ",s[1:].translate())

print("3rd frame: ",s[2:].translate())

rev=s.reverse_complement()

print("4th frame: ",rev.translate())

print("5th frame: ",rev[1:].translate())

print("6th frame: ",rev[2:].translate())