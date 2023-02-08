l1=[12,34,563456,345,45,"CATGACTG",-9.76]

print(l1)

print(type(l1))

#access the values from list

print(l1[1])

print(l1[-1])

#range of values from list

print(l1[2:5])

print(l1[1:])

print(l1[:7])

print(l1[:])

print(l1[::-1])

#length of list

print(len(l1))

#change valye ofrom list

l1[1]=1000

print(l1)

#add element in list

l1.append("gc")

print(l1)

#add element at speciifc location
l1.insert(0,1999)

print(l1)

#remove

l1.remove(45)
print(l1)

l1.pop(1)
print(l1)

#loop through list

for i in l1:
    print(i)

#join two list

l2=[12,34,45,45,99]
l3=["acc","SFF","RRT"]

l=l2+l3
print(l)


#sort

l4=[12.23,34,45,1,45,56,12,67,45]

l4.sort()

print(l4)

l4.sort(reverse=True)

print(l4)

print(max(l4))

print(min(l4))

print(sum(l4))

l4.remove(34)
print(l4)

print(l4)
l4.remove(12)
print(l4)


print(l4)
d1=dict(l4)
print(d1)
