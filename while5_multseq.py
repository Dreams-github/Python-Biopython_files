import collections

s1="""ATTGGACTTTGAACTGATAAAAACGGAAGCCAAACATCCTGCCACTTTAAGGAAGTATTGTATAGAGGCA
AAGCTGACCAACACAACTACAGCATCTCGTTGCCCAACACAAGGAGAACCCAGCCTAAATGAAGAACAGG
ACAAAAGGTTTGTCTGCAAACACTCCATGGTAGACAGAGGATGGGGAAATGGATGCGGATTGTTTGGAAA
GGGAGGCATCGTGACCTGTGCAATGTTCACATGCAAAAAGAACATGGAAGGAAAAGTCGTGCAACCAGAA
AACTTGGAGTATACCATTGTGATAACACCTCACTCAGGGGAAGAGAATGCAGTCGGAAATGACACAGGAA
AACACGGCACGGAAATTAAAGTAACGCCACAGAGTTCCATCACAGAAGCGGAACTGACAGGCTATGGCAC
TGTCACGATGGAATGCTCTCCGAGAACGGGCCTCGACTTTAATGAGATGGTGTTGCTGCAAATGGAAGAC
AAGGCTTGGCTGGTGCACAGGCAATGGTTCTTAGACCTGCCGTTACCATGGCTGCCCGGAGCAGACAAAC
AAGGATCAAATTGGATACAGAAGGAGACATTGGTCACTTTCAAAAATCCCCATGCGAAGAAACAGGATGT
TGTTGTTTTAGGATCCCAAGAAGGGGCCATGCATACAGCACTCACAGGGGCCACGGAAATCCAGATGTCA
TCAGGAAACTTACTG"""
print("combination of patterns are CTTT,AAAA,AAC,GCCC")

p="Y"

while(p=="Y"):
    comb=str(input("enter the pattern from the above combination"))
    if(comb=="CTTT"):
        print(comb,"present")
        c=s1.count(comb)
        print("count of CTTT",comb,c)
    elif(comb=="AAAA"):
        print(comb,"present")
        c=s1.count(comb)
        print("count of AAAA", comb,c)
    elif(comb=="AAC"):
        print(comb,"present")
        c=s1.count(comb)
        print("count of AAC",comb,c)
    elif(comb=="GCCC"):
        print(comb, "present")
        c = s1.count(comb)
        print("count of GCCC", comb, c)
    else:
        print("the entered pattern is not form the above combination")
        print(comb, " not present")
        c = s1.count(comb)
        print("count of diff comb", comb, c)
    p = str(input("Do you want to continue?? (Y/N) ")).upper()

