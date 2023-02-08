from Bio.Seq import Seq

S1=Seq("ATTGCAGTGACCCCGTAGGTAAAGTGACCGATCGCTGGATCGGCGTAGGCGAAATGCTGGGGGGGGGGGACCCAGGATGATCGTGG")

print(S1)

prot=S1.translate()

print(prot)

l1=[]
for i in range(len(prot)):
    if(prot.startswith("*",i)):
        #print(i)
        l1.append(i)

print(l1)

print(prot[:l1[2]])