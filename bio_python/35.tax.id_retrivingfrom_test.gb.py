from Bio import SeqIO

data=SeqIO.parse("test.gb","genbank")

for i in data:
    for id in i.features:
        if(id.type=="source"):
            id="".join(id.qualifiers.get("db_xref"))
            print(i.id,id)

