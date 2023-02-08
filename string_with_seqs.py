s1="TAGTCGATGCTGCTGTCGATCTGTAGTCGTGATGCTGACTGCTGCTGTGCTGTAGTCGAGTCG"
s2="ATCGTGCTCGTAGCTCGTAGCTGCTGCTGCTGCTAGCGATGTCGTCGCGATGTGCT"
s3="ACTGCTCGTCGTCGTCGTCGATGCGTGCTGCTCTGTGTCGTCG"

s=s1+s2+s3
print("added seq",s)

#access the char from seq
print(s[2])

print(s[44],s[49])

print(s[-1])

#access the range of chars

print(s[3:7])

print(s[5:100])

print(s[6:])

print(s[:55])

print(s[:])

#length of seq

print(len(s))

#count A T G C

print("count A",s.count("A"))

print("count T",s.count("T"))

print("count G",s.count("G"))

print("count C",s.count("C"))

#GC %

g=s.count("G")
c=s.count("C")

pGC=((g+c)/len(s))*100
print("GC %",pGC)

a=s.count("A")
t=s.count("T")

pAT=((a+t)/len(s))*100
print("AT %",pAT)

#dna to rna

print(s.replace("T","U"))

#split

ss1="AGTCGATCCG TCAGTCGCGCT ATCGTCGTACG ATCTCGCT"

print(ss1.split(" "))

#remove spaces

ss2="   TACTCGTCGCGT   TACCGATGTGCTGC   ATGCTGATCG   TACGTCATGCTGC   TACGTGCGT"

print(ss2.replace("   ",""))

#strip

ss3="  TCATCGTCGCAGCGTAGTATCTGA   "

print(ss3.strip())

#FIND

print(s.find("TCG"))

print(s.count("TCG"))

#startswith

print(s.startswith("TAG"))

#endswith

print(s.endswith("GC"))