#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "par_sim.h"
#include "int_lennardjones.h"


double fintx = 0;
double finty = 0;
double fintxpair = 0;
double fintypair = 0;

/**function for the calculation of the intra-particle force arising
 * from a Lennard-Jones Potential. The potential is described by two
 * parameters, position of minimum LJ_MIN and depth of minimum  LJ_EPS.
 * The potential is truncated at INT_CUTOFF. 
 */
/*inline double intforce_lj(double dist1d, double dist2d){

	return(LJPREFAC*(LJMINPOW*pow(dist2d,-14)-pow(dist2d,-8))*dist1d);

}

inline double intforce(double dist1d, double dist2d){

	return intforce_lj(dist1d, dist2d);

}*/

void specs_int(double f_cut){
	
	FILE *outpspecs;
	outpspecs = fopen("muovert_specs.dat", "a");
	fprintf(outpspecs, "\n\nParameters of Lennard_Jones particle-particle interaction:\n\nRadius of Particles' Hardcore: %.4lf\nPotential cut-off length: %.2lf\nValue of int-force at cut-off length: %.5lf\nDepth of potential minimum epsilon: %.2lf\nPosition of potential minimum: %.2lf\nPower of 6 of minimum's position: %.4e\n\n\n", R_INT, INT_CUTOFF, f_cut, EPS_L, LJMIN, LJMINPOW);

	fclose (outpspecs);



}


void copycode_int(){
char copycode[200];

  sprintf(copycode, "cp ../int_lennardjones.* ./");
  system(copycode);

}

char* prfx_int(){
char *a = "LJ_pot";
return a;
}
