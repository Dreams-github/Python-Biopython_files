f1=open("file1.fasta","r")
f2=open("file2.fasta","r")

data1=f1.read()
data2=f2.read()

if(data1==data2):
    print("both files having same content")
else:
    print("not same")