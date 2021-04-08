from math import *
import os

#static force along x-axis
f = 0.64
#sample number of interacting particles
sn = 5
#number of tasks if parallelized with MPI
numtasks = 2

#content of run_script (specific for albeniz cluster)
str_run_script =\
"#!/bin/bash\n\
##$ -N sn_{0}_{1}\n\
##$ -o $JOB_NAME-$JOB_ID.log\n\
##$ -j y\n\
#$ -cwd\n\
#$ -V\n\
# echo 'Running on' `hostname`\n\
mpirun -n  {2} main_brownconf {3:.3f} {4}\n\
scp your_output_files.dat  nodo00:$SGE_O_WORKDIR"

#ask if mpicc for parallelization with mpi is installed
mpi_valid = os.popen("which mpicc").read()

#if shell output is empty, mpi is not installed
if mpi_valid != '':
	print(mpi_valid)
	mpi_flag = 0
else:
	mpi_flag = 1


#while loop over range of forces to initialize
#jobs with respective forcefind = f
find = f
while fabs(find) <=  fabs(f):
#while fabs(find) <=  2:
  scriptname = "LJ_sn_{0}_F_{1:.2f}".format(sn,find) 
  file_name = open(scriptname, "w")
  file_name.write(str_run_script.format(f, sn, numtasks, find, sn))
  file_name.close()

  print("f={:.3f}".format(find))
  print("setnumb={}".format(sn))
  runscript = 'qsub ' + scriptname
  #runscript = 'qsub -pe mpi 1 ' + scriptname
  os.system(runscript)  
  find = find*2


