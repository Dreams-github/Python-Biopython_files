from Bio.Seq import Seq

s1=Seq("""CGACTCCAAGGAGAAACTGGAGAGAGCATGGACTAGCAGCTTGACAGCACACCTAAAGCTCTAGAACAAA
AAGAAGCAAATACACTCAAGAAGAGTAGAAGTTAAGAAATAATCAAACTCAGAGTTGGAATCGACCAAGT
AGAAACAAAAAGGACTATACAAAGAATCAACAAAAGCAGGAGAGGTTCTTTGAGAAAATCAACAAGATAG
ATAAATCCTTAGCCAGACTAACCAGAGGGCACAGAAAGTGTATTCAAATTAACAAAACCAGAAATGAAAA
GGGAGACATAACAAAAGAATCTGAGCAAAAAAAAAAAA
""")
#print(s1)
#print(len(s1))
s2=s1.replace("\n","")
print(s2)
print(len(s2))
print(s2.translate())
print(s2.translate(stop_symbol="@"))
print(s2.translate(to_stop=True))
print(s2.translate(to_stop=False))
#print(s2.translate(cds=True))
print(s2.translate(table=1))
print(s2.translate(table=2))
print(s2.translate(table="Yeast Mitochondrial",cds=True))
print(s2.translate(gap="-"))