import os

if(os.path.exists("mm.txt")):
    os.remove("mm.txt")
else:
    print("file is not present")