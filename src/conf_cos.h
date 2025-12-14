/**ensure that only one header for spatial
 * confinement is included in code*/
#ifndef HEADER_CONF
#define HEADER_CONF

#include <math.h>
#include <stdbool.h>
#include "sim_config.h"

#define AMP (1.0/(2.0*M_PI))
#define MAX_HALF_WIDTH (2*AMP+BOTTLENECK_WIDTH)

#define CHECKN 5
#define K_COS (2.0*M_PI/L_CONF)

#define XMAX (L_CONF/4.0-R_CONF)
#define SQRT_SHIFT R_CONF*sqrt(1+AMP*AMP*K_COS*K_COS)

/**
 *********************************************************
 *
 * External types and variables
 *
 *********************************************************
 */


/**
 *********************************************************
 *
 * Internal functions
 *
 *********************************************************
 */

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
static inline double yuef_cos(double x, double y){

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
		return (BOTTLENECK_WIDTH + AMP + AMP*sin(K_COS*x));
	}
	/*If particle position is well below confinement let particle pass*/
	if((BOTTLENECK_WIDTH + AMP + AMP*sin(K_COS*x) - fabs(y)) > SQRT_SHIFT){
		return(2*fabs(y));
	}
	/*If ypos is larger than ymax don't let particle pass and leave function*/
	ymax = BOTTLENECK_WIDTH+AMP+AMP*sin(K_COS*XMAX)-R_CONF*sqrt(1+AMP*AMP*K_COS*K_COS*cos(K_COS*XMAX)*cos(K_COS*XMAX));	
	if(fabs(y) >= ymax){
		return(ymax);
	}

	/*if particle is in the vicinity of boundaries, start to scan the particle's surface*/ 
	xcirc = R_CONF;
	PosValid = true;
	do{
		xcirc = xcirc - 2*R_CONF/CHECKN;
		xt = xcirc + x;
		yu = BOTTLENECK_WIDTH + AMP + AMP*sin(K_COS*xt);
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

/**
 *********************************************************
 *
 * External functions
 *
 *********************************************************
 */

/*function handler that provides generic interface to main*/
static inline double CONF_yuef(double x, double y){
	return yuef_cos(x, y);
}

/*extern double yuef_ext(double x, double y);
extern double yuef_cos(double x, double y);*/
extern void CONF_specs(char *file_confparams);
extern void CONF_copycode();
extern char *CONF_prfx();

#endif
