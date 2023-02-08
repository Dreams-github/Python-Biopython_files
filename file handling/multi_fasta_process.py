import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
f=open("ls_orchid.fasta","r")
f12=open("out(1).csv", "w")
data=f.read()

print("id"+"\t"+"countA"+"\t"+"countT"+"\t"+"countG"+"\t"+"countC"+"\t"+"GC%"+"\t"+"AT%",file=f12)

#print(data)

sp_data=data.split(">")
print(sp_data)
cntA=[]
cntT=[]
cnG=[]
cnC=[]
GC=[]
AT=[]
header=[]
seqs=[]
for s in sp_data:
    if(s!=""):
        #print(s)
        head=s[:s.find("\n")]
        print(head)
        seq=s[s.find("\n"):]
        #print(seq)
        nseq=seq.replace("\n","")
        print(nseq)
        cA=nseq.count("A")
        cT=nseq.count("T")
        cG=nseq.count("G")
        cC=nseq.count("C")
        print("count of A",cA)
        print("count of T",cT)
        print("count of G",cG)
        print("count of c",cC)
        PGC=((cG+cC)/len(nseq))*100
        PAT=((cA+cT)/len(nseq))*100
        print("GC %",PGC)
        print("AT %",PAT)

        id = head[:head.find(" ")]
        header.append(id)
        seqs.append(nseq)
        cntA.append(cA)
        cntT.append(cT)
        cnG.append(cG)
        cnC.append(cC)
        GC.append(PGC)
        AT.append(PAT)

        out=id+"\t"+str(cA)+"\t"+str(cT)+"\t"+str(cG)+"\t"+str(cC)+"\t"+str(PGC)+"\t"+str(PAT)
        print(out,file=f12)


df=pd.DataFrame(
    {
        "head":header,
        "seqs":seq,
        "countA":cntA,
        "countT":cntT,
        "countG":cnG,
        "countC":cnC,
        "GC%":GC,
        "AT%":AT

    }
)

print(df)

writer=ExcelWriter("pandas_excel.xlsx")
df.to_excel(writer,"my_data",index=False)
writer.save()