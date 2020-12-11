/**
 * Module to handle a confinement with a cosine shape 
 * (see e.g. Phd Thesis Steffen Martens 2013 page 90).
 * Function with boundary for septate channel 
 * is in header file for inlining and thus a
 * better computational performance
 */
#include <stdio.h>
#include <stdlib.h>
#include "par_sim.h"
#include "conf_cos.h"
#include "results_transport.h"


/**function to write confinement specific information in specs file*/
void CONF_specs(){

    FILE *outpspecs;
    outpspecs = fopen("muovert_specs.dat", "a");		
    fprintf(outpspecs, "\n\nParameters of 'Cosine'-Confinement:\n\nChannel Length L: %.1lf\nBottleneck Half-Width B: %.2lf\nAmplitude A: %.2lf\nChannel's max. Half-Width: %.2lf\nChannel's wave number: %.2lf\nParticle Radius: %.10lf\nBottrad B/R: %.3lf\n\nBinwidth 1d x-Histogram: %lf\nBinwidht 1d y-Histogram: %lf\nBinwidth 2d Histogram: %lf\n\n", L, 
																	 B, 
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
