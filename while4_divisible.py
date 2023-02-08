#100 to 900
cD=0
a=100
while(a<=900):
    if(a%3==0):
        print("the value divisible by 3",a)
        cD=cD+1
    else:
        print("not",a)
    a=a+1
print("the count of divisible",cD)