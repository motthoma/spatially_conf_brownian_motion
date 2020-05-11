import os
import sys


SERVER =  "motthoma@161.116.80.211"

def read_out_content_server():
	#function that provides list of all zip files on 
	#server that contain code_push_key in name
	str_ls_server = "ssh " + SERVER +  " ls" 
	ls_out_server = os.popen(str_ls_server).read()
	list_out_server = ls_out_server.split('\n')
	list_code_zips = []
	for item in list_out_server:
		if ('.zip' in item):
 			#remove .zip ending to easily compare file names in following
			list_code_zips.append(item.split(".zip")[0])
	return list_code_zips

def print_list_index(print_list, message):
    #print list elements with corresponding index
    #message string is printed before list elements    
    print(message)
    print('ind\t name')
    for ind, item in enumerate(print_list):
    	print('{}:\t {}.zip'.format(ind, item))
    print('\n')

def get_name_zip(zip_selection, action = 'download'):
    #get name of zip file to be downloaded from or uploaded to server
    args = list(sys.argv)
    if len(args) > 1:
    	zip_in = args[1]
    else:
    	file_select_request = "choose the zip file to be " + action + "ed by typing name or index: "
        #use try except construction to enable python 2 and 3 compatibility
    	try:
    		zip_in = raw_input(file_select_request) 
    	except:
        	zip_in = input(file_select_request)
  
    #check if input is index number:
    if zip_in[0].isdigit():
    	zip_selected_suf = zip_selection[int(zip_in)]
    else:
    	zip_selected_suf = zip_in

   
    if zip_selected_suf.endswith(".zip"):
    	zip_selected = zip_selected_suf.split(".zip")[0]
    else:
    	zip_selected = zip_selected_suf

    return zip_selected
 
def get_selected_zip(zip_selected, list_code_zips):
    #function that validates name of selected zip and picks out
    #suitable zip from server

    #truncate name of wanted zip to length of longest available zip name
    #skip .zip suffix to avoid a partial truncation of .zip suffix 
    max_len_zip = max([len(item) for item in list_code_zips])
    if len(zip_selected) > max_len_zip:
    	zip_selected = zip_selected[:max_len_zip]

    #check if wanted zip is in list of available ones
    chosen_zips = []
    for zip_file in list_code_zips:
    	if zip_selected == zip_file:
        	chosen_zips.append(zip_selected)
   
    #if none of the available zips matches,
    #find the ones where name fits partially
    if chosen_zips == []:
	    chosen_zips = []
	    for zip_file in list_code_zips:
	    	if zip_selected in zip_file:
                   	chosen_zips.append(zip_file)
    return chosen_zips

def check_zip_list(list_zips):
    #stop if name was not clear so that several zips could fit
    if len(list_zips) == 0:
    	raise Exception('no matching zip file was found to input. Check printed list of available zips')
    elif len(list_zips) > 1:
    	raise Exception('too many matching zip files were found. Check printed list of available zips')

    chosen_file = list_zips[0]
    if chosen_file.endswith(".zip"):
    	return chosen_file
    else:
    	return chosen_file + ".zip"


def download_selected_zip(zip_file_name):
    #function that downloads selected zip from server
    string_copy = "scp -r " + SERVER + "://home/motthoma/" + zip_file_name + ' ./'
    os.system(string_copy)

def unzip_file(zip_file_name):
    string_unzip = "unzip " + zip_file_name + " -d ./" + zip_file_name.split(".zip")[0]
    os.system(string_unzip)

if __name__ == "__main__":
  
    list_code_zips = read_out_content_server() 
#    list_code_zips = ['code_brownian_2.zip', 'code_brownian_int_17_04_2020_23_25_19.zip', 'code_brownian_int_from_albeniz.zip', 'code_brownian_int.zip', 'code_brownian_test.zip']   
    
    print_list_index(list_code_zips, '\nlist of found .zip files on ' + SERVER + ':\n')

    zip_selected = get_name_zip(list_code_zips)
    print("\ndemanded zip file: {}\n".format(zip_selected))
 
    list_chosen_zips = get_selected_zip(zip_selected, list_code_zips)
    print("suitable zips found on {}: {}\n".format(SERVER, list_chosen_zips))
    
    chosen_zip = check_zip_list(list_chosen_zips)
    print("zip file selected for download: {}\n".format(chosen_zip))  
  
    download_selected_zip(chosen_zip)	
    print("{} successfully downloaded\n".format(chosen_zip))

    unzip_file(chosen_zip) 
