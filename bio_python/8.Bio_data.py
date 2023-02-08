from Bio.Data import CodonTable

st_table=CodonTable.unambiguous_dna_by_id[1]

print(st_table)

mt_table=CodonTable.unambiguous_dna_by_id[2]
print(mt_table)

print(type(CodonTable.unambiguous_dna_by_id))

print(CodonTable.unambiguous_dna_by_name.keys())


mm=CodonTable.unambiguous_dna_by_name["Mold Mitochondrial"]
print(mm)

print(mm.start_codons)

print(mm.stop_codons)