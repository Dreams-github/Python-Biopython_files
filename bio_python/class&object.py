class seq:
    head=""
    ss=""

    def __init__(self,h,s):
        self.head=h
        self.ss=s


    def countAT(self):
        a=self.ss.count("A")
        t=self.ss.count("T")
        print("count of AT",a+t)

    def countGC(self):
        g=self.ss.count("G")
        c=self.ss.count("C")
        print("count of GC",g+c)

obj1=seq(">h","ATGCGTGTCGTAAGTGCCCAAGTGGGGTGG")

print(obj1.ss)

obj1.countAT()
obj1.countGC()

obj2=seq(">h2","CGTGAGTGCGTAGTGGAAGTCGTAAAGTGAC")

print(obj2.head,obj2.ss)

obj2.countGC()
obj2.countAT()