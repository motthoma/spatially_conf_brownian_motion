/* Module that contains functions to handle a confinement given
 * by a septated channel according to: 
 * Borromeo, Marchesoni Chem. Phys. 2010, doi = "https://doi.org/10.1016/j.chemphys.2010.03.022"
 *
 * Function with boundary for septate channel 
 * is in header file for inlining and thus a
 * better computational performance
*/ 
#include <stdio.h>
#include <stdlib.h>
#include "code_handling.h"
#include "conf_sept.h"
#include "results_transport.h"

void CONF_specs(){
    
    char fname_intparams [60]; 
    sprintf(fname_intparams,"parameters_confinement.dat");

    DestPaths.fname_confparams = malloc(256);    
    snprintf(DestPaths.fname_confparams,
             800,
             "%s/%s",
             DestPaths.fullpath, fname_intparams);

    FILE *outpspecs;
    outpspecs = fopen(DestPaths.fname_confparams, "a");		
    fprintf(outpspecs,
            "\n\nParameters of Confinement 'Septated Channel':\n\n"
		    "Channel Length L_CONF: %.1lf\n"
		    "Bottleneck Half-Width B: %.2lf\n"
		    "Amplitude AMP: %.2lf\n"
		    "Channel's max. Half-Width: %.2lf\n"
		    "Particle Radius: %.2lf\n"
		    "Bottrad B/R: %.3lf\n\n"
		    "Binwidth 1d x-Histogram: %lf\n"
		    "Binwidht 1d y-Histogram: %lf\n"
		    "Binwidth 2d Histogram: %lf\n\n", 
             L_CONF, 
             BOTTLENECK_WIDTH, 
             AMP, 
             MAX_HALF_WIDTH, 
             R_CONF, 
             BOTTRAD, 
             histparams.binx, 
             histparams.biny, 
             histparams.bin2d);
    
    
fclose(outpspecs);
}

void CONF_copycode(){
    CODEHAND_copy_file_to_dest("conf_sept.c");
    CODEHAND_copy_file_to_dest("conf_sept.h");
}

char* CONF_prfx(){
    char *a = "sept";
    return a;
}
