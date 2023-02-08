seq="ACGTATGCTGAGTGACGTGATGGATGCGTAGTGGACCGTGTAGTGACGTCGCCC"

start_location=0

for i in range(len(seq)):
    if(seq[i:i+3]=="ATG"):
        start_location=i
        break
print(start_location)

for j in range(start_location,len(seq),3):
    codon=seq[j:j+3]
    print(codon)