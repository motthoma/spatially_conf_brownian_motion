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
#include "par_sim.h"
#include "conf_sept.h"
#include "results_transport.h"

void CONF_specs(char *file_confparams){

    FILE *outpspecs;
    outpspecs = fopen(file_confparams, "a");		
    fprintf(outpspecs, "\n\nParameters of Confinement 'Septated Channel':\n\n"
		    "Channel Length L: %.1lf\n"
		    "Bottleneck Half-Width B: %.2lf\n"
		    "Amplitude AMP: %.2lf\n"
		    "Channel's max. Half-Width: %.2lf\n"
		    "Particle Radius: %.2lf\n"
		    "Bottrad B/R: %.3lf\n\n"
		    "Binwidth 1d x-Histogram: %lf\n"
		    "Binwidht 1d y-Histogram: %lf\n"
		    "Binwidth 2d Histogram: %lf\n\n", 
		    				     L, 
						     B, 
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
    char copycode[200];

     sprintf(copycode, "cp ../conf_sept.* ./");
     system(copycode);

}

char* CONF_prfx(){
    char *a = "sept";
return a;
}
