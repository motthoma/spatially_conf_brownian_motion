"""

Script that helps to create makefile
for the compilation of brownian motion
code.

It presents the possible files for 
confinements and interactions and composes
the makefile based on the user's choices.

"""
import os
import dtool_pull_code_from_server as pull_code


def get_module_files(prefix, message):
	#function that creates and prints list of c-modules of respective type
	files = [f for f in os.listdir('.') if os.path.isfile(f) and f.endswith('.c')]
	conf_files = [f for f in files if f.startswith(prefix)]

	pull_code.print_list_index(conf_files, message)	
	
	return conf_files


def get_compiler():
	"""
	function that contains list of suggested compilers and takes user
	input to select compiler
	"""	
	list_compilers = ['gcc', 'clang', 'icc', 'mpicc']
	message = "\n\nList of suggested compilers:\n"
	pull_code.print_list_index(list_compilers, message)	
	
	message = "\nchoose compiler used for compilation if in list. If not, type enter:"
	compiler = select_file(list_compilers, message)
	if compiler == []:
		compiler = pull_code.scan_input("\ntype in wanted compiler:")	
	return compiler
	
def select_file(list_select, scan_message):	
        #function that takes user input to select c-module from list	
	
	select = pull_code.scan_input(message)
	
	if select.isdigit():
		file_sel = [list_select[int(select)]]
	else:
		file_sel = pull_code.get_selected_file(select, list_select)
	
	return file_sel[0]

def adapt_line(line, new_part):
	#function that cuts part right of '=' away and replaces it by new_part
	left_part = line.split('=')[0]
	new_line = left_part + '= ' + new_part + '\n'	
	return new_line

def adapt_makefile(compiler, conf_file, int_file):
	
	makefile_temp = open('makefile_temp', 'w')
	makefile_old =  open('makefile', 'r')

	for line in makefile_old:
		if line.startswith('CC'):
			new_line = adapt_line(line, compiler)
			makefile_temp.write(new_line)
		
		elif line.startswith('CONF'):
			new_line = adapt_line(line, conf_file)
			makefile_temp.write(new_line)
		
		elif line.startswith('INT'):
			new_line = adapt_line(line, int_file)
			makefile_temp.write(new_line)
		else:
			makefile_temp.write(line)

	makefile_temp.close()
	makefile_old.close()
	os.system('cp makefile_temp makefile')


def write_conf_int_header(conf_file, int_file):
#function that creates header file that contains
#includes header for confinement and inter-particle
#interactions
	conf_string = conf_file.split('.c')[0]
	conf_string += '.h'
	int_string = int_file.split('.c')[0]
	int_string += '.h'

	out_h = open("comp_gen_header.h", 'w')

	out_h.write('\n')

	header_interact = int_string
	out_h.write('#include "{}"\n'.format(header_interact))

	header_conf = conf_string
	out_h.write('#include "{}"\n'.format(header_conf))

	out_h.close()

def call_make_file():
	os.system("make")


if __name__ == "__main__":

	message = "\n\nList of possible confinments:\n"
	conf_list_files = get_module_files('conf', message)
	message = "\nchoose module for confinement to be linked in compilation:"
	conf_file = select_file(conf_list_files, message)
	print(conf_file)

	message = "\n\nList of possible interactions:\n"
	int_list_files = get_module_files('int', message)		
	message = "\nchoose module for interaction to be linked in compilation:"
	int_file = select_file(int_list_files, message)
	print(int_file)	

	compiler = get_compiler()
	print(compiler)	
	
	adapt_makefile(compiler, conf_file, int_file)
	
	write_conf_int_header(conf_file, int_file)

	call_make_file()
