import os

if __name__ == "__main__":
    
    zip_file_name = "zip code_brownian_int.zip"
    string_zip = zip_file_name + " *.c *.h *.py makefile"
    os.system(string_zip)

    string_copy = "scp -r " + zip_file_name + " motthoma@161.116.80.211://home/motthoma" 
    os.system(string_copy)

