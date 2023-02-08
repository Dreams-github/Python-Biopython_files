s="TAGCTGCAGCTGCTAGTGGAGACGGCGCATGCTGACTTGACTGTAGCTGACTGACTGCATGCATCATGCTGCTGCGCATGCTGAACCA"

print("original seq",s)

#ATGC
s=s.replace("A","a")  #aTGC
s=s.replace("T","A")  #aAGC
s=s.replace("a","T")  #TAGC
s=s.replace("G","g")  #TAgC
s=s.replace("C","G")  #TAgG
s=s.replace("g","C")  #TACG

print("complement of seq",s)

seq="TGCATGCTGGACTGATGACTTCGTACGTTACGTACGGACCCTGATAATGACTGACTTGAC"

print("original sequence ;",seq)

seq=seq.replace("A","a")
seq=seq.replace("T","A")
seq=seq.replace("a","T")
seq=seq.replace("G","g")
seq=seq.replace("C","G")
seq=seq.replace("g","C")

print("cpmlimentary sequence:",seq)



