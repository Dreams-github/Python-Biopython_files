from Bio.Seq import Seq
from Bio import SeqUtils

sdata=Seq("""CGACTCCAAGGAGAAACTGGAGAGAGCATGGACTAGCAGCTTGACAGCACACCTAAAGCTCTAGAACAAA
AAGAAGCAAATACACTCAAGAAGAGTAGAAGTTAAGAAATAATCAAACTCAGAGTTGGAATCGACCAAGT
AGAAACAAAAAGGACTATACAAAGAATCAACAAAAGCAGGAGAGGTTCTTTGAGAAAATCAACAAGATAG
ATAAATCCTTAGCCAGACTAACCAGAGGGCACAGAAAGTGTATTCAAATTAACAAAACCAGAAATGAAAA
GGGAGACATAACAAAAGAATCTGAGCAAAAAAAAAAAA""")

print(len(sdata))

s1=sdata.replace("\n","")
print(s1)
print(len(s1))

print(s1.count("A"))
print(s1.count("T"))
print(s1.count("G"))
print(s1.count("C"))
print("GC %",SeqUtils.GC(s1))
AT=100-SeqUtils.GC(s1)
print("AT %",AT)#print("molecular weight of s1",SeqUtils.molecular_weight(s1))
print(s1.back_transcribe())
print("six frames",SeqUtils.six_frame_translations(s1))
#b=s1.back_transcribe()
#print("six frames",SeqUtils.six_frame_translations(b))