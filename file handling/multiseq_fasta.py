f=open("sequence (4).fasta","r")
f11=open("out.csv","w")
data=f.read()
#print(data)

s_data=data.split(">")
#print(s_data)

print("id"+"\t"+"countA"+"\t"+"countT"+"\t"+"countG"+"\t"+"countC"+"\t"+"GC%"+"\t"+"AT%",file=f11)
for s in s_data:
    if(s!=""):
        #print(s)
        head=s[:s.find("\n")]
        print(head)
        seq=s[s.find("\n"):]
        #print(seq)
        s1=seq.replace("\n","")
        print(s1)
        cA=s1.count("A")
        cT=s1.count("T")
        cG=s1.count("G")
        cC=s1.count("C")
        print("count of A",cA)
        print("count of T",cT)
        print("count of G",cG)
        print("count of C",cC)
        gc=cG+cC
        at=cA+cT
        PGC=gc/len(s1)*100
        PAT=at/len(s1)*100
        print("GC %",PGC)
        print("AT %",PAT)

        id=head[:head.find(" ")]

        out=id+"\t"+str(cA)+"\t"+str(cT)+"\t"+str(cG)+"\t"+str(cC)+"\t"+str(PGC)+"\t"+str(PAT)
        print(out,file=f11)

