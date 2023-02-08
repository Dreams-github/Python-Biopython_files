s1={34,78,56,90,"AR",-99,8.7}

print(s1)

print(type(s1))

print(len(s1))

for i in s1:
    print(i)

print(s1)

s1.add(89)
print(s1)

s1.remove(-99)
print(s1)

s2={23,45,"ADTFG",56,-156}

s=s1.union(s2)
print(s)


l1=[89,156,354,456,456,78,98,89]
print(l1)

s4=set(l1)
print(s4)