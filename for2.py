
flag="Y"
while(flag=="Y"):
    m=int(input("enter the number"))

    for a in range(1,11):
        r=m*a
        print(m,"x",a,"=",r)
    flag=str(input("Do you want to continue Y/N? ")).upper()
