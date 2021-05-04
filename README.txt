********************************************************************
*								   *
*								   *
*								   *
*	Project code to simulate interacting and non-interacting   *
*	brownian particles in 2-dimensional periodic channels      *
*								   *
*								   *
*								   *
********************************************************************


OVERVIEW OF PROJECT DIRECTORIES


This is a description of the simulation code for interacting and- non interacting particles. The code can be used to simulate 
non-interacting and interacting particles in periodic 2-dimensional channels. The particles perform an overdamped brownian motion,
which is simulated with a stochastic Euler procedure. The particles can be exposed to external static forces and do not 
experience long-range hydrodynamic interactions. Periodic boundary conditions are employed while the particle density is 
kept constant if interacting particles are simulated.



The project contains three directories:

1. code -> contains the code, makefile and latest simulation results. Each simulation run starts with the creation of a directory where the code and data is stored.

2. doxygen -> contains the doxygen code analysis. type doxygen to run the Doxyfile as well as the code analysis in form of html files and a pdf created by pdf latex.

3. profiling -> contains an environment to profile the runtime of the code with gprof. Open the README therein for furthe information. 

DETAILS ON COMPILATION AND CODE EXECUTION

In order to allow simulations of several models, the code contains different modules accounting for various channel shapes (prefix 'conf_') and modules for hard-core and lennard-jones interactions (prefix 'int_').

If MPI is installed (such as on Albeniz cluster), the code can also be prallelized.

The choice of the specific model is made by linking the wanted modules when compiling the code. i
The code can be compiled with a Makefile, while dtool_create_makefile.py is a script that guides the user during the choices to be made and creates a Makefile. During the process, also the wanted compiler can be chosen.

The code can be executed by the execution of 
	
		./main_brownconf {force} {ensemble_size} {mpi_jobs}

with the external static force, the number of interacting particles (not the number of all simulated particles!) and the number of parallel mpi jobs. If non-interacting particles are to be simulated, simply choose any type of interaction and set the ensemble size to one.


FURTHER TOOLS

The directocry 'code' contains further tools:

1. masterinteract.py -> script that can be used to initiate several jobs on albeniz
2. dtool_push/pull_code_from/to_server.py -> scripts that can be used to zip the code and push/pull it to/from albeniz
