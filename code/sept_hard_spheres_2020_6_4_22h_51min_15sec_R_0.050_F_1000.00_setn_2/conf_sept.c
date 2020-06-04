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


void specs_conf(double binx, double biny, double bin2d){

    FILE *outpspecs;
    outpspecs = fopen("muovert_specs.dat", "a");		
    fprintf(outpspecs, "\n\nParameters of Confinement 'Septated Channel':\n\nChannel Length L: %.1lf\nBottleneck Half-Width B: %.2lf\nAmplitude AMP: %.2lf\nChannel's max. Half-Width: %.2lf\nParticle Radius: %.2lf\nBottrad B/R: %.3lf\n\nBinwidth 1d x-Histogram: %lf\nBinwidht 1d y-Histogram: %lf\nBinwidth 2d Histogram: %lf\n\n", L, B, AMP, MAX_HALF_WIDTH, R_CONF, BOTTRAD, binx, biny, bin2d);
    
    
fclose(outpspecs);
}

void copycode_conf(){
    char copycode[200];

     sprintf(copycode, "cp ../conf_sept.* ./");
     system(copycode);

}

char* prfx_conf(){
    char *a = "sept";
return a;
}
