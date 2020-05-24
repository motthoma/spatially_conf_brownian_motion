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


def get_compiler(mpicc):
	"""
	function that contains list of suggested compilers and takes user
	input to select compiler
	"""	
	list_compilers = ['gcc', 'clang', 'icc']
	if mpicc == True:
		list_compilers.append('mpicc')
	message = "\n\nList of suggested compilers:\n"
	pull_code.print_list_index(list_compilers, message)	
	
	message = "\nchoose compiler by index or name used for compilation if in list. If not, type enter:"
	compiler = select_file(list_compilers, message)
	if compiler == []:
		compiler = pull_code.scan_input("\ntype in wanted compiler:")	
	return compiler
	
def select_file(list_select, scan_message):	
        #function that takes user input to select c-module from list	
	
	select = pull_code.scan_input(scan_message)
	
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
	#adaption of make file in order to link the respecitve
	#module for the chosen confinement and inter-particle
	#interaction which have been included in header file
	#linked to main function of c-Cocde
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
	os.system('rm makefile_temp')


def write_conf_int_header(conf_file, int_file, mpi_flag):
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
	
	if mpi_flag == True:
		out_h.write('#include "mpi.h"\n')
		out_h.write('#define MPI_ON\n')
	out_h.close()

def call_make_file():
	os.system("make")

def check_for_mpicc():
	#ask if mpicc for parallelization with mpi is installed
	mpi_valid = os.popen("which mpicc").read()

	#if shell output is empty, mpi is not installed
	if mpi_valid != '':
		mpi_flag = True
	else:
		mpi_flag = False

	return mpi_flag

if __name__ == "__main__":

	list_title_conf = "\n\nList of possible confinments:\n"
	conf_list_files = get_module_files('conf', list_title_conf)
	request_conf = "\nchoose module for confinement to be linked in compilation (ind or name):"
	conf_file = select_file(conf_list_files, request_conf)
	print(conf_file)

	list_title_int = "\n\nList of possible interactions:\n"
	int_list_files = get_module_files('int', list_title_int)		
	request_int = "\nchoose module for interaction to be linked in compilation (ind or name):"
	int_file = select_file(int_list_files, request_int)
	print(int_file)	

	mpicc_flag = check_for_mpicc()

	compiler = get_compiler(mpicc_flag)
	print(compiler)	
	
	adapt_makefile(compiler, conf_file, int_file)
	
	write_conf_int_header(conf_file, int_file, mpicc_flag)

	call_make_file()
