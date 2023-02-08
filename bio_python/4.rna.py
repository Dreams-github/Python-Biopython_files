from Bio.Seq import Seq

s1=Seq("AUCGUAAGUUUUCGUAGCUAGCUGAUCGGUUCGAUCGAAGUCGGAAUUUU")
print(s1)

print(s1.complement())
print(s1.back_transcribe())
print(s1.translate())