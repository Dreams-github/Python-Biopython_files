class sequence:
    header=""
    seq=""

    def __init__(self,H,S):
        self.header=H
        self.seq=S

    def countGC(self):
        g=self.seq.count("G")
        c=self.seq.count("C")
        print("count GC",g+c)

    def countAT(self):
        a = self.seq.count("A")
        t = self.seq.count("T")
        print("count AT", a + t)

    def rna(self):
        print("rna is",self.seq.replace("T","U"))

ob1=sequence(">h1","TGACTGTAGTCGTACTGTCGTGCATGTCCGCTGTGCTTCGTGCTGCTGCTCTGCTGTCGTC")

ob1.countAT()

ob1.rna()

ob2=sequence(">h222","ATCTGTTAGTCGTAGCTGCTGATGTCGTAGCGTACGTGCAGCTGCATCTCGATGCTAGCTCGTAGCTGTAG")

ob2.countGC()