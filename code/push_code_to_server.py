import os
from datetime import datetime
import pull_code_from_server as pull_code


SERVER = pull_code.SERVER

def replace_warning(zip_name, list_zips):
    #show warning if user is about to overwrite existing file on server
    if zip_name.endswith(".zip"):
        zip_name = zip_name.split(".zip")[0]
    else:
        zip_name = zip_name

    if zip_name in list_zips:
	print("Warning!\n A file with the same name as the file to be uploadet already exists at destination!\n")

	over_write_warning = "Do you want to overwrite existing file? Type [y/n]:"
	#use try except construction to enable python 2 and 3 compatibility
        try:
                answer = raw_input(over_write_warning)
        except:
                answer = input(over_write_warning)

	if 'y' in answer:
		pass
        else:
	        raise Exception("push process is stopped")

def create_zip(zip_file_name):
    #create zip file
    string_zip = zip_file_name + " *.c *.h *.py makefile"
    os.system(string_zip)

def push_zip(zip_file_name):
    #copy zip file to SERVER

    string_copy = "scp -r " + zip_file_name + SERVER + "://home/motthoma" 
    os.system(string_copy)

if __name__ == "__main__":
    
    list_code_zips = pull_code.read_out_content_server()   
 
    pull_code.print_list_index(list_code_zips, '\nlist of found .zip files on ' + SERVER + ':\n') 

    zip_name = pull_code.get_name_zip(list_code_zips, 'upload')

    replace_warning(zip_name, list_code_zips)

    if zip_name.endswith(".zip"):
       zip_name_suf = zip_name
    else:
       zip_name_suf = zip_name + ".zip"	

    create_zip(zip_name_suf)
 
    push_zip(zip_name_suf)
 
