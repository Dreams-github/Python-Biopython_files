t1=(12,34,45,34,-43.34,"ACTCGC",-99)

print(t1)

print(type(t1))

#access the elements

print(t1[1])

print(t1[2])

print(t1[3:7])

print(t1[:4])

print(t1[::-1])

print(len(t1))

for i in t1:
    print(i)

t2=(5,8,19,24)
t3=(25,54,45,42)
t4=t2+t3
print(t4)

l=list(t4)
print(l)

l.append(67)
print(l)

l.insert(5,52)
print(l)

t=tuple(l)
print(t)

sum(t)
print(sum(t))
print(max(t))
print(min(t))
print(t)


l1=list(t)
print(l1)

l1.sort()
print(l1)

l1.sort(reverse=True)
print(l1)