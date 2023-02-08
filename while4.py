s="ATGGACGGTTATGGACGGTACAACATGGACGGTAATGGGCAGTAACGATGGGCTATGGTAATGGGTAACGGCATGGGCAGTAATGGAGACATGGGTAAGAC"
print("combinations are ATGG,GCA,GTA")

flag="Y"

while(flag=="Y"):
    nuc=str(input("enter the pattern from above combination. "))
    if(nuc=="ATGG"):
        print(nuc,"present")
        c=s.count(nuc)
        print("count of",nuc,c)
    elif(nuc=="GCA"):
        print(nuc,"present")
        c=s.count(nuc)
        print("count of",nuc,c)
    elif(nuc=="GTA"):
        print(nuc,"present")
        c=s.count(nuc)
        print("count of",nuc,c)
    else:
        print("entered pattern is not from above combination")
        print(nuc,"present")
        c=s.count(nuc)
        print("count of",nuc,c)

    flag=str(input("Do you want to continue?? (Y/N) "))