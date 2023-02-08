s1="welcome"
s2="to"
s3="pune"

print(s1)
print(s2)
print(s3)

#concatenate the strings

s=s1+" "+s2+" "+s3
print(s)

#welcome to pune
#0123456789......  index

#access the characters

print(s[1])

print(s[5],s[8])

print(s[-1])

#range of characters

print(s[0:9])

print(s[5:10])

print(s[5:])

print(s[:10])

print(s[:])

print(s[::-1])

#length of strings

print(len(s))

#count char

print(s.count("e"))

#replace

print(s.replace("e","E"))

#upper

print(s.upper())

#lower

print(s.lower())

print(s.capitalize())

#split

print(s.split("e"))

#strip

ss="  python  "
print(ss)

print(ss.strip())

#lstrip

print(ss.lstrip())

#r strip

print(ss.rstrip())

#find

print(s.find("e"))

#startswith

print(s.startswith("wel"))

#endswith

print(s.endswith("ne"))

#islower

print(s.islower())

#isupper

print(s.isupper())

#format

p1="this is {} from {}".format("shubham","sangli")

print(p1)