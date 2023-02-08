f=open("blast_result.txt","r")

for s in f:
    #print(s)
    sp=s.strip().split("\t")
    print(sp)
    length=float(sp[3])
    mismatch=float(sp[4])
    gap=float(sp[5])
    EML=length-mismatch-gap
    QC=((EML)/float(sp[-1]))*100
    print("QC",QC)