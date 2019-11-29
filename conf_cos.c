#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "par_sim.h"
#include "conf_cos.h"



/**function to account for a confinement with cosine-shape.
 * An analytical value of the effecitve confinement at the respective particle 
 * position is only available for small spherical particles whose boundary curvature
 * is smaller than the curvature of the confinement (see e.g. Phd Thesis
 * Steffen Martens 2013 page 90).
 * Therefore, if the particle at position (x,y) is in the vicinity 
 * of the boundary, the positions of the particle's surface are scanned and checked
 * if they lie beyond the boundary.
 * If particle is within confinement, a virtual value for the effective confinement
 * yuef of yuef = 2*y is provided to ensure that yuef > y. 
 * Otherwise, yuef = 0.5*y is provided to ensure that y < yuef.
 */
inline double yuef_cos(double x, double y){

/**relative position within the particle, where the y-value of the particle
 * surface is checked if particle is within the confinement*/
double xcirc;
/**absolute position of scan*/
double xt;
/**value of confinement*/ 
double yu;
/**height of particle's surface at scan position xt*/
double ycirc;
/**maximal value of confinement*/
double ymax;
/**counter to break up routine if particle is outside of confinement*/
bool PosValid;
	
	/*If no finite sized particles are considered simply analytical 
         * formula is provided.
         */ 
	if(R_CONF == 0){
		return (B + AMP + AMP*sin(K_COS*x));
	}
	/*If particle position is well below confinement let particle pass*/
	if((B + AMP + AMP*sin(K_COS*x) - fabs(y)) > SQRT_SHIFT){
		return(2*fabs(y));
	}
	/*If ypos is larger than ymax don't let particle pass and leave function*/
	ymax = B+AMP+AMP*sin(K_COS*XMAX)-R_CONF*sqrt(1+AMP*AMP*K_COS*K_COS*cos(K_COS*XMAX)*cos(K_COS*XMAX));	
	if(fabs(y) >= ymax){
		return(ymax);
	}

	/*if particle is in the vicinity of boundaries, start to scan the particle's surface*/ 
	xcirc = R_CONF;
	PosValid = true;
	do{
		xcirc = xcirc - 2*R_CONF/CHECKN;
		xt = xcirc + x;
		yu = B + AMP + AMP*sin(K_COS*xt);
		ycirc = sqrt(R_CONF*R_CONF-xcirc*xcirc)+fabs(y);
		/*leave routine if particle is outside of confinement and deliver yeff 
                 *thats below y-value of particle*/ 
		if (ycirc >= yu){
			PosValid = false;
			return(0.5*fabs(y));
		}
	}while((PosValid == false) && (xcirc > -R_CONF));
	/*if particle is within confinement, deliver yeff thats twice of particle coordinate y*/
	return(2*fabs(y));
}

/*function handler that provides generic interface to main*/
inline double yuef_ext(double x, double y){
	
	return yuef_cos(x, y);
}




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
