import os

if __name__ == "__main__":
    
    zip_file_name = "zip code_brownian_2.zip"
    string_zip = zip_file_name + " *.c *.h *.py makefile"
    os.system(string_zip)
    string_copy = "scp -r " + "motthoma@161.116.80.211://home/motthoma/" + zip_file_name ' ./'
    os.system(string_copy)

