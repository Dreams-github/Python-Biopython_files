s="ATGGCCGCGATGGTTGGTATGGCGCGATATGGCGCGTGGATGGCGCGTTTATGGCCCGCGTTTATGGCCATGGCCGTT"

pattern="ATGG"

if(pattern in s):
    print(pattern,"present")
    c=s.count(pattern)
    print("count of",pattern,c)
else:
    print(pattern,"not presebt")