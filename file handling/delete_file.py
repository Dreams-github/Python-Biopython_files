import os

if(os.path.exists("new_file.py")):
    os.remove("new_file.py")
else:
    print("new_file.py")