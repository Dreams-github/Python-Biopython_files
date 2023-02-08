from Bio.Seq import Seq
s1=Seq("""CGGGTTGTAGCTGGCCTTTAACGAGGTATGTGCACGCCTGGCTCATCCACTCTCAACCTCTGTGCACTTT
ATGTAAGAAACGGTGTAAGCCAGCTATTCATTAGTTGGTAATAAGCCTTTCTTATGTTTACTACAAACGC
TTCAGTTATAGAATGTTTACTGTGTATAACACAATTATATACAACTTTCAGCAACGGATCTCTTGGCTCT
CGCATCGATGAAGAACGCAGCGAAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATC
TTTGAACGCACCTTGCACTCCTTGGTATTCCGAGGAGTATGCCTGTTTGAGTCTCATGGAATTCTCAACC
CCTAAATTTTGTAATGAAGTTTATTGGGCTTGGACTTGTAGGTTGTGTCGGCTTCTATTCAACTCCTCTG
AAAA""")

print("original seq",s1)
print(len(s1))
s2=s1.replace("\n","")
print("single line seq",s2)
print(len(s2))
print("complement seq",s2.complement())
print("reverse complement seq",s2.reverse_complement())
print("seq transcribe",s2.transcribe())
print("seq translate",s2.translate())