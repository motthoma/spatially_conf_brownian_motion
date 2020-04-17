import os
from datetime import datetime

if __name__ == "__main__":

    now = datetime.now()
    str_now = now.strftime("%d_%m_%Y_%H_%M_%S")

    zip_file_name = "zip code_brownian_int_{}.zip".format(str_now)
    string_zip = zip_file_name + " *.c *.h *.py makefile"
    os.system(string_zip)

    string_copy = "scp -r " + zip_file_name + " motthoma@161.116.80.211://home/motthoma" 
    os.system(string_copy)

