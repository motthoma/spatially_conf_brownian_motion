import os

if __name__ == "__main__":

    zip_file_name = "code_brownian_int.zip"
    string_copy = "scp -r " + "motthoma@161.116.80.211://home/motthoma/" + zip_file_name + ' ./'
    os.system(string_copy)

    string_unzip = "unzip " + zip_file_name + " -d ./code_unzip"
    os.system(string_unzip)
