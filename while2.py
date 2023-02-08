#40 to 234  even odd
cE=0
cD=0
a=40
while(a<=243):
    if(a%2==0):
        print("the even values",a)
        cE=cE+1
    else:
        print("the odd value",a)
        cD=cD+1
    a=a+1
print("the count of even",cE)
print("the count of odd",cD)



