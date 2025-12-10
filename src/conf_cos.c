/**
 * Module to handle a confinement with a cosine shape 
 * (see e.g. Phd Thesis Steffen Martens 2013 page 90).
 * Function with boundary for septate channel 
 * is in header file for inlining and thus a
 * better computational performance
 */
#include <stdio.h>
#include <stdlib.h>
#include "conf_cos.h"
#include "results_transport.h"


/**function to write confinement specific information in specs file*/
void CONF_specs(char *file_confparams){

    FILE *outpspecs;
    outpspecs = fopen(file_confparams, "a");		
    fprintf(outpspecs, "\n\nParameters of 'Cosine'-Confinement:\n\n"
		     	"Channel Length L_CONF: %.1lf\n"
			"Bottleneck Half-Width B: %.2lf\n"
			"Amplitude A: %.2lf\n"
			"Channel's max. Half-Width: %.2lf\n"
			"Channel's wave number: %.2lf\n"
			"Particle Radius: %.10lf\n"
			"Bottrad B/R: %.3lf\n\n"
			"Binwidth 1d x-Histogram: %lf\n"
			"Binwidht 1d y-Histogram: %lf\n"
			"Binwidth 2d Histogram: %lf\n\n", 
							L_CONF, 
							BOTTLENECK_WIDTH, 
						 	AMP, 
							MAX_HALF_WIDTH,
							K_COS, 
							R_CONF, 
							BOTTRAD, 
							histparams.binx, 
							histparams.biny, 
							histparams.bin2d);
			    
    
fclose(outpspecs);
}

void CONF_copycode(){
char copycode[200];

  sprintf(copycode, "cp ../conf_cos.* ./");
  system(copycode);

}

/**function to provide prefix in name of working directory*/
char* CONF_prfx(){
char *a = "cos";
return a;
}
