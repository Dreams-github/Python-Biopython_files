d1={"geneID":"NC899","SEQ":"ATGGCCGTTT","length":7887}
print(d1)
print(type(d1))

#access the value by key
print(d1["geneID"])

#add new k-v pair

d1["GC%"]=60

print(d1)

d1["AT%"]=75
print(d1)

#length
print(len(d1))

#loop thr

for j in d1:
    print(j,d1[j])

for a,b in d1.items():
    print(a,b)

#keys

print(d1.keys())
print(d1.values())

#copy

d2=d1.copy()
print(d2)

