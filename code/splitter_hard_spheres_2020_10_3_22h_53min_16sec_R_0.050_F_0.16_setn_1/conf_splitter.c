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
#include "par_sim.h"
#include "conf_splitter.h"



void CONF_specs(double binx, double biny, double bin2d){

    FILE *outpspecs;
    outpspecs = fopen("muovert_specs.dat", "a");		
    fprintf(outpspecs, "\n\nParameters of 'Splitter'-Confinement:\n\nChannel Length L: %.1lf\nBottleneck Half-Width B: %.2lf\nSlope M: %.2lf\nChannel's max. Half-Width: %.2lf\nParticle Radius: %.2lf\nBottrad B/R: %.3lf\nLp: %.2lf\n\nBinwidth 1d x-Histogram: %lf\nBinwidht 1d y-Histogram: %lf\nBinwidth 2d Histogram: %lf\n\n", L, B, M, MAX_HALF_WIDTH, R_CONF, BOTTRAD, Lp, binx, biny, bin2d);
    
    
fclose(outpspecs);
}

void CONF_copycode(){
char copycode[200];

  sprintf(copycode, "cp ../conf_splitter.* ./");
  system(copycode);

}

/**
 * function that provides name of confinement
 */
char* CONF_prfx(){
	char *a = "splitter";
	return a;
}
