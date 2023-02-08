s1= "today"
s2= "is"
s3= "monday"
print(s1)
print(s2)
print(s3)

s=s1+" "+s2+" "+s3
print(s)

print(len(s))

print(s.count("y"))

print(s.upper())

print(s.lower())



print(s.capitalize())

print(s.split("d"))
print(s)

print(s.strip())

b=" arraygen  "

print(b)

print(b.strip())

print(b.lstrip())
print(b.rstrip())
print(s)

print(s.startswith("mon"))
print(s.isupper())
#find

print(s.find("y"))

print(s[:].find("y"))

p="my {} is {} from today".format("session","starting")
print(p)

print(p[5])
print(p[:10])

print(p[5:])

#string with sequence

seq1= "ATGCTGACTGACTGACTTGTCTATGGGTCGTAG"
seq2="ATGCGTGCTGAGTCGTAGCCCGTAGTVGTATG"
seq3="ATGCATGTCGTAGCGTGTAGGTCAGTGCGGGCTGAACGTGATCTTTGGACCCTGA"

seq=seq1+seq2+seq3
print("adding sequences :",seq)

print(seq.find("GAC"))

print(len(seq))

print("count GC",seq.count("GC"))
#count of ATGC

print("count of A :",seq.count("A"))
print("count of T :",seq.count("T"))
print("count of G :",seq.count("G"))
print("count of C :",seq.count("C"))

#GC%

g= seq.count("G")
c= seq.count("C")

PGC=((g+c)/len(seq))*100
print("GC%=",PGC)

#AT%

PAT=((seq.count("A")+seq.count("T"))/len(seq))*100
print("AT%=",PAT)

ss="""GTCACGCTAGAAGAGCCGTCTGCTGACTGCACGTGTGTGTGCACACTCGTGTGCATGGCTGTGAACTGGA
ATGTGTGACTGTGACCTTATGGCTGCCGCACGCCTCTGCCTCTCCCTGCTGCTCCTGTCCACCTGCGTGG
CTCTGTTACTACAGCCACTGCTGGGTGCCCAGGGAGCCCCACTGGAGCCAGTGTACCCAGGGGACAATGC
CACACCAGAGCAGATGGCCCAGTATGCAGCTGATCTCCGTAGATACATCAACATGCTGACCAGGCCTAGG
TATGGGAAAAGACACAAAGAGGACACGCTGGCCTTCTCGGAGTGGGGGTCCCCGCATGCTGCTGTCCCCA
GGGAGCTCAGCCCGCTGGACTTATAATGCCACCTTCTGTCTCCTACGACTCCATGAGCAGCGCCAGCCCA
GCTCTCCCCTCTGCACCCTTGGCTCTGGCCAAAGCTTGCTCCCTGCTCCCACACAGGCTCAATAAAGCAA
GTCAAAGCCA"""

print(ss)

print(len(ss))

new=ss.replace("\n","")
print(new)

print(len(new))

print(ss.replace("T","U"))

print(ss.lower())
print(ss.isupper())
print(ss.startswith("gtc"))

print(ss.endswith("GCCA"))





