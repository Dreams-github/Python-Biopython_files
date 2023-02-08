#100 to 500
s=0
b=0
i=100
while(i<=500):
    if(i%2==0):
        print("even",i)
        s=s+i
    else:
        print("odd",i)
        b=b+i
    i=i+1
print("the sum of even",s)
print("the sum of odd",b)

