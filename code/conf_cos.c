#include <stdio.h>
#include <stdlib.h>
#include "par_sim.h"
#include "conf_cos.h"



/**function to write confinement specific information in specs file*/
void specs_conf(double binx, double biny, double bin2d){

    FILE *outpspecs;
    outpspecs = fopen("muovert_specs.dat", "a");		
    fprintf(outpspecs, "\n\nParameters of 'Cosine'-Confinement:\n\nChannel Length L: %.1lf\nBottleneck Half-Width B: %.2lf\nAmplitude A: %.2lf\nChannel's max. Half-Width: %.2lf\nChannel's wave number: %.2lf\nParticle Radius: %.10lf\nBottrad B/R: %.3lf\n\nBinwidth 1d x-Histogram: %lf\nBinwidht 1d y-Histogram: %lf\nBinwidth 2d Histogram: %lf\n\n", L, B, AMP, MAX_HALF_WIDTH, K_COS, R_CONF, BOTTRAD, binx, biny, bin2d);
    
    
fclose(outpspecs);
}

void copycode_conf(){
char copycode[200];

  sprintf(copycode, "cp ../conf_cos.* ./");
  system(copycode);

}

/**function to provide prefix in name of working directory*/
char* prfx_conf(){
char *a = "cos";
return a;
}
