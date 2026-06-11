#include <stdio.h>
#include <stdlib.h>
#include "int_lennardjones.h"
#include "code_handling.h""


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

void INT_specs(){
	
    char fname_intparams [60]; 
    snprintf(fname_intparams,
             sizeof fname_intparams,
             "parameters_particle_interaction.dat");
    DestPaths.fname_intparams = malloc(256);    
    snprintf(DestPaths.fname_intparams,
             256,
             "%s/%s",
             DestPaths.fullpath, fname_intparams);

	FILE *outpspecs;
	outpspecs = fopen(DestPaths.fname_intparams, "a");
	fprintf(outpspecs,
            "\n\nParameters of Lennard_Jones particle-particle interaction:\n\n"
            "Radius of Particles' Hardcore: %.4lf\n"
            "Potential cut-off length: %.2lf\n"
            "Value of int-force shift value: %.5lf\n"
            "Depth of potential minimum epsilon: %.2lf\n"
            "Position of potential minimum: %.2lf\n"
            "Power of 6 of minimum's position: %.4e\n\n\n",
            R_INT, 
            INT_CUTOFF, 
            FSHIFT, 
            EPS_L, 
            LJMIN, 
            LJMINPOW);

	fclose (outpspecs);
}

void INT_copycode(){
    CODEHAND_copy_file_to_dest("int_lennardjones.c");
    CODEHAND_copy_file_to_dest("int_lennardjones.h");
}

char * INT_prfx(){
	char *a = "LJ_pot";
    return a;
}
