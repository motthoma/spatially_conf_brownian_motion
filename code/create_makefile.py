"""

Script that helps to create makefile
for the compilation of brownian motion
code.

It presents the possible files for 
confinements and interactions and composes
the makefile based on the user's choices.

"""
import os
import pull_code_from_server as pull_code


def get_confinement_files():
	#function that creates and prints list of confinment c-modules
	files = [f for f in os.listdir('.') if os.path.isfile(f) and f.endswith('.c')]
	conf_files = [f for f in files if f.startswith('conf')]

	print("\n\nList of possible confinments:\n")	
	for ind, element in enumerate(conf_files):
		print("{}: {}".format(ind, element))
	
	
	return conf_files

def get_interaction_files():
	#function that creates and prints list of interaction c-modules

	files = [f for f in os.listdir('.') if os.path.isfile(f) and f.endswith('.c')]
	int_files = [f for f in files if f.startswith('int')]
	
	print("\n\nList of possible interactions:\n")	
	for ind, element in enumerate(int_files):
		print("{}: {}".format(ind, element))
	
	return int_files

def select_file(select, list_select):	
	
	if select.isdigit():
		file_sel = [list_select[int(select)]]
	else:
		file_sel = pull_code.get_selected_file(select, list_select)
	
	return file_sel		

if __name__ == "__main__":

	conf_list_files = get_confinement_files()
	conf_select = pull_code.scan_input("\nchoose file for confinement to be linked in compilation:")
	conf_file = select_file(conf_select, conf_list_files)
	print(conf_file)

	int_list_files = get_interaction_files()		
	int_select = pull_code.scan_input("\nchoose file for interaction to be linked in compilation:")
	int_file = select_file(int_select, int_list_files)
	print(int_file)	
