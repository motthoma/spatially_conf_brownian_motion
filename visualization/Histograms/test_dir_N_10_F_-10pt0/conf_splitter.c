/* Module that contains functions to handle a confinement given
 * by a saw tooth like profile used for the entropic splitter 
 * (see Motz et al. J. Chem. Phys. 2014).
 *
 * Function with boundary for septate channel 
 * is in header file for inlining and thus a
 * better computational performance
*/ 
#include <stdio.h>
#include <stdlib.h>
#include "code_handling.h"
#include "conf_splitter.h"
#include "results_transport.h"


void CONF_specs(){
    char fname_intparams [60]; 
    snprintf(fname_intparams,
             sizeof fname_intparams,
             "parameters_confinement.dat");

    DestPaths.fname_confparams = malloc(256);    
    snprintf(DestPaths.fname_confparams,
             800,
             "%s/%s",
             DestPaths.fullpath, fname_intparams);

    FILE *outpspecs;
    outpspecs = fopen(DestPaths.fname_confparams, "a");		
    fprintf(outpspecs, "\n\nParameters of 'Splitter'-Confinement:\n\n"
		    "Channel Length L_CONF: %.1lf\n"
			"Bottleneck Half-Width BOTTLENECK_WIDTH: %.2lf\n"
			"Slope M: %.2lf\n"
			"Channel's max. Half-Width: %.2lf\n"
			"Particle Radius: %.2lf\n"
			"Bottrad B/R: %.3lf\n"
			"Lp: %.2lf\n\n"
			"Binwidth 1d x-Histogram: %lf\n"
			"Binwidht 1d y-Histogram: %lf\n"
			"Binwidth 2d Histogram: %lf\n\n", 
            L_CONF, 
            BOTTLENECK_WIDTH, 
            SLOPE, 
            MAX_HALF_WIDTH, 
            R_CONF, 
            BOTTRAD, 
            Lp, 
            histparams.binx, 
            histparams.biny, 
            histparams.bin2d);

    
fclose(outpspecs);
}

void CONF_copycode(){
    CODEHAND_copy_file_to_dest("conf_splitter.c");
    CODEHAND_copy_file_to_dest("conf_splitter.h");
}

/**
 * function that provides name of confinement
 */
char* CONF_prfx(){
	char *a = "splitter";
	return a;
}
