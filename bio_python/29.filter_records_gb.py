from Bio import SeqIO

data=SeqIO.parse("test.gb","genbank")

for s in data:
    c=0
    for i in s.features:
        if(i.qualifiers.get("collection_date")):
            #print(s.id,i.qualifiers.get("collection_date"))
            date="".join(i.qualifiers.get("collection_date"))
            if(date!="None"):
                print(s.id,date)
                c=True
        if(i.type=="CDS" and c==True):
            s1 = i.location._start.position
            s2 = i.location._end.position
            print(s1,s2)
            print(s.seq[s1:s2])


