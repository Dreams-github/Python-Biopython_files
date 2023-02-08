def path_fun(path):
    f=open(path,"r")

    cA=0
    cT=0
    cG=0
    cC=0

    for s in f:
        if(">" in s):
            print(s)
        else:
            print(s)
            cA=cA+s.count("A")
            cT=cT+s.count("T")
            cG=cG+s.count("G")
            cC=cC+s.count("C")
    print("count of A",cA)
    print("count of T",cT)
    print("count of G",cG)
    print("count of C",cC)

path_fun("C:/Users/avksa/PycharmProjects/pythonProject/file handling/sequence (2).fasta")
path_fun(r"C:\Users\avksa\OneDrive\Desktop\sequence (5).fasta")