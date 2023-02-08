from Bio.Seq import Seq

s=Seq("ATTGCAGTGACCCCGTAGGTAAAGTGACCGATCGCTGGATCGGCGTAGGCGAAATGCTGGGGGGGGGGGACCCAGGATGATCGTGG")

print(s)

print(s.complement())
print(s.reverse_complement())
print(s.transcribe())
#print(s.translate())

if(len(s)%3==0):
    print("the seq is paired: ",s.translate())
elif(len(s)%3==1):
    print("the seq remain one bp: ",s[:-1].translate())
elif(len(s)%3==2):
    print("the seq remain two bp: ",s[:-2].translate())

