#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "par_sim.h"
#include "conf_sept.h"


/**function to provide effective boundary of confinement with the shape of 
 * septated channels (see e.g. Marchesoni J. Chem. Phys. 2010).
 * Effective boundary can only be represented piece wise: In the vicinity of the
 * bottleneck, yueff is a circle with radius R_CONF. In the other sections, we
 * have shifted straight lines parallel to the originial boundary.
 */
inline double yuef_sept(double x, double y){

/*If particle is in central cylinder of width 2B,
 * an evaluation of the eff. boundary is not necessary.
 * Provide value of confinement that ensures y < yueff.
 */
if(fabs(y) < B - R_CONF){
	return(MAX_HALF_WIDTH);
}

/*evaluate effective confinement piece wise. Confinement
 * is given by a circle of Radius R_CONF in the vicinity 
 * of the bottlenecks.*/
if ((R_CONF <= x) && (x < L - R_CONF)){
	return (B + AMP - R_CONF);
}
	
        if (x <= R_CONF){
		return (B - sqrt(R_CONF_SQ - x*x));
	}

		if ((L-R_CONF <= x) && (x <=L)) {
			return (B - sqrt(R_CONF_SQ - (x - L)*(x - L)));
		}

return(10);
}

/*function handler that provides generic interface to main*/
inline double yuef_ext(double x, double y){
	
	return yuef_sept(x, y);
}

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
