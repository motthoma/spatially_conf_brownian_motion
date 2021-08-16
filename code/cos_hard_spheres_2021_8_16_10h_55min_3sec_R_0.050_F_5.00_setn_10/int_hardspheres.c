#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "int_hardspheres.h"


double fintx = 0;
double finty = 0;
double fintxpair = 0;
double fintypair = 0;

/*inline int intforce(double disti, double dist){

	return(0);

}*/

void INT_specs(char *intspecs){
	
	double f_cut;
  	f_cut = INT_force(INT_CUTOFF, INT_CUTOFF);
	
	FILE *outpspecs;
	outpspecs = fopen(intspecs, "a");
	fprintf(outpspecs, "\n\nParameters of hard-sphere particle-particle interaction:\n\n"
			   "Radius of Particles' Hardcore: %.2lf\n"
			   "Potential cut-off length: %d\n"
			   "Value of int-force at cut-off length: %.5lf\n\n", R_INT, 
			   						      INT_CUTOFF, 
									      f_cut);

	fclose (outpspecs);



}

void INT_copycode(){
char copycode[200];

  sprintf(copycode, "cp ../int_hardspheres.* ./");
  system(copycode);

}

char* INT_prfx(){
	char *a = "hard_spheres";
	return a;
}
