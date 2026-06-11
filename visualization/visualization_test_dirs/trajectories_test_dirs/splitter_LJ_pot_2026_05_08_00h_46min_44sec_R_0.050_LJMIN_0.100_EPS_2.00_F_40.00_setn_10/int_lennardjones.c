#include <stdio.h>
#include <stdlib.h>
#include "int_lennardjones.h"
#include "code_handling.h"

double fintx = 0;
double finty = 0;
double fintxpair = 0;
double fintypair = 0;

/** Writes specifications of the interaction to the specification output file. */
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

/**
 * copies the interaction specific code to the destination folder 
 * for documentation purposes
 */
void INT_copycode(){
    CODEHAND_copy_file_to_dest("int_lennardjones.c");
    CODEHAND_copy_file_to_dest("int_lennardjones.h");
}

/**
 * provides the interaction specific prefix in name of working directory
 */
char * INT_prfx(){
	char *a = "LJ_pot";
    return a;
}
