my_rna="GCAUAUGUUCAUAGACGAGUUAAUCAGCUGA"

star_location=0

for i in range(len(my_rna)):
    if(my_rna[i:i+3]=="AUG"):
        star_location=i
        break


print(star_location)

for j in range(star_location,len(my_rna),3):
    codon=my_rna[j:j+3]
    print(codon)