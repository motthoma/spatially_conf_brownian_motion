#include <stdio.h>
#include <stdlib.h>
#include "code_handling.h"
#include "int_hardspheres.h"


double fintx = 0;
double finty = 0;
double fintxpair = 0;
double fintypair = 0;

/*inline int intforce(double disti, double dist){

	return(0);

}*/

void INT_specs(){
	
	double f_cut;
  	f_cut = INT_force(INT_CUTOFF, INT_CUTOFF);
	
    char fname_intparams [60]; 
    sprintf(fname_intparams, "parameters_particle_interaction.dat");
    DestPaths.fname_intparams = malloc(256);    
    snprintf(DestPaths.fname_intparams,
             800,
             "%s/%s",
             DestPaths.fullpath, fname_intparams);

	FILE *outpspecs;
	outpspecs = fopen(DestPaths.fname_intparams, "a");
	fprintf(outpspecs, "\n\nParameters of hard-sphere particle-particle interaction:\n\n"
			   "Radius of Particles' Hardcore: %.2lf\n"
			   "Potential cut-off length: %d\n"
			   "Value of int-force at cut-off length: %.5lf\n\n", R_INT, 
			   						      INT_CUTOFF, 
									      f_cut);

	fclose (outpspecs);
}

void INT_copycode(){
    CODEHAND_copy_file_to_dest("int_hardspheres.c");
    CODEHAND_copy_file_to_dest("int_hardspheres.h");
}

char* INT_prfx(){
	char *a = "hard_spheres";
	return a;
}
