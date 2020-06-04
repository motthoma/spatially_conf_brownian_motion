"""
script that profiles code of interacting brownian particles project. 
The code is copied from 'code' directory to local one, compiled, executed
and profiled
"""

import os
import sys

#adjustable parameters
#choose name of executable, type of confinement and inter-particle interaction
executable = "main_brownconf"
confinement = "sept"
interaction = "hardspheres"

#choose value of external force, size of interacting samples
#and number of prallel threads
ext_force = "1"
setnumb = "1"
numb_mpi_tasks = "0"

if __name__ == "__main__":
	#copy source code to local directory
	code_path = "../code"
	copy_str_cfiles = "cp " + code_path + "/*.c ./" 
	copy_str_header = "cp " + code_path + "/*.h ./"
	copy_str_masterfile = "cp " + code_path + "/masterinteract.py ./" 
	print(copy_str_cfiles)
	os.system(copy_str_cfiles)
	os.system(copy_str_header)
	os.system(copy_str_masterfile)

	#compile source code with debugging -p -pg flags for profiling
	compile_str = "gcc -Wall main_brownconf.c par_sim.c conf_{0}.c int_{1}.c -lgsl -lgslcblas -lm  -o {2} -p -pg".format(confinement, interaction, executable)
	print(compile_str)
	os.system(compile_str)

	#execute source code to create gmon.out file
	execute_str = "./{0} {1} {2} {3}".format(executable, ext_force, setnumb, numb_mpi_tasks)
	os.system(execute_str) 

	#perform profiling
	#check if gmon.out file exists
	if os.path.isfile('gmon.out'):
	    pass
	else:
	    raise Exception("gmon.out does not exist!")

	output = "gprof_result_conf_{0}_int_{1}_force_{2}_setnumb_{3}.txt".format(confinement, interaction, ext_force, setnumb)
	profile_str = "gprof ./{0} gmon.out > {1}".format(executable, output)
	os.system(profile_str)
