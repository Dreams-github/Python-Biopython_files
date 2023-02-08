from Bio.Seq import Seq

s1=Seq("TACGCTGCATCTTTTGCATGCATCATGCAGTCGATCGCTAGTGCTGCTATGATGTTCGTACGTTATCAGAGCTGTCTCTGCTGGACTGCGGGCACCTAG")

print("orginal seq",s1)

#complement

print("complement of seq",s1.complement())

#reverse commplement

print("reverse complement",s1.reverse_complement())

#transcribe

print("transcribe",s1.transcribe())

#translate

print("translate ",s1.translate())