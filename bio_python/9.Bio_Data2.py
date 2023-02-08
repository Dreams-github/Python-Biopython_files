from Bio.Data import CodonTable

tb1=CodonTable.unambiguous_dna_by_id[1]

print(tb1)
print(tb1.start_codons)
print(tb1.stop_codons)

print(type(CodonTable.ambiguous_dna_by_id))
print(CodonTable.ambiguous_dna_by_id.keys())
print(CodonTable.ambiguous_dna_by_name.keys())

tb2=CodonTable.unambiguous_dna_by_name["Yeast Mitochondrial"]
print(tb2)



