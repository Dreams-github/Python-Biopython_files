from Bio.Seq import Seq

ss="""CGGGTTGTAGCTGGCCTTTAACGAGGTATGTGCACGCCTGGCTCATCCACTCTCAACCTCTGTGCACTTT
ATGTAAGAAACGGTGTAAGCCAGCTATTCATTAGTTGGTAATAAGCCTTTCTTATGTTTACTACAAACGC
TTCAGTTATAGAATGTTTACTGTGTATAACACAATTATATACAACTTTCAGCAACGGATCTCTTGGCTCT
CGCATCGATGAAGAACGCAGCGAAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATC
TTTGAACGCACCTTGCACTCCTTGGTATTCCGAGGAGTATGCCTGTTTGAGTCTCATGGAATTCTCAACC
CCTAAATTTTGTAATGAAGTTTATTGGGCTTGGACTTGTAGGTTGTGTCGGCTTCTATTCAACTCCTCTG
AAAA"""

print(ss)

s1=ss.replace("\n","")
print(s1)
s2=Seq(s1)

print("frame 1: ",s2.translate())
print("frame 2: ",s2[1:].translate())
print("frame 3: ",s2[2:].translate())

rev=s2.reverse_complement()

print("frame 4: ",rev.translate())
print("frame 5: ",rev[1:].translate())
print("frame 6: ",rev[2:].translate())