seq="""AGTTCTAGTTCTATCTTTTTCCTGCATCCTGTCTGGAAGTTAGAAGGAAACAGACCACAGACCTGGTCCC
CAAAAGAAATGGAGGCAATAGGTTTTGAGGGGCATGAGGACGGGGTTCAGCCTCCAGGGTCCTACACACA
AATCAGTCAGTGGCCCAGAAGACCCCCCTCGGAATCGGAGCAGGGAGGATGGGGAGTGTGAGGGGTATCC
TTGATGCTTGTGTGTCCCCAACTTTCCAAATCCCCGCCCCCGCGATGGAGAAGAAACCGAGACAGAAGGT
GCAGGGCCCACTACCGCTTCCTCCAGATGAGCTCATGGG"""
print(len(seq))
s1=seq.replace("\n","")
print(s1)
print(len(s1))



c=0
f="Y"
while(f=="Y"):
    pattern = str(input("enter the pattern")).upper()
    for i in range(len(s1)):
        if(s1.startswith(pattern,i)):
            print("starting location",pattern,i,"ending location",pattern,i+len(pattern)-1)
            c=c+1
    print("count of the pattern",c)

    f=str(input("do you want to continue?? Y/N")).upper()

