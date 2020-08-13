/**
 * ensure that only one header for spatial
 * confinement is included in code
 * Function with boundary for septate channel 
 * is in header file for inlining and thus a
 * better computational performance
 */
#ifndef HEADER_CONF
#define HEADER_CONF

#include <math.h>
#include "par_sim.h"

/*define parameters which are not defined in
  par_sim.h header file
*/  
#define AMP 1.0*L
#define MAX_HALF_WIDTH (AMP+B)
#define R_CONF_SQ R_CONF*R_CONF

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

/**function to provide effective boundary of confinement with the shape of 
 * septated channels (see e.g. Marchesoni J. Chem. Phys. 2010).
 * Effective boundary can only be represented piece wise: In the vicinity of the
 * bottleneck, yueff is a circle with radius R_CONF. In the other sections, we
 * have shifted straight lines parallel to the originial boundary.
 */
static inline double yuef_sept(double x, double y){

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

/**
 *********************************************************
 *
 * External functions
 *
 *********************************************************
 */

/*function handler that provides generic interface to main*/
static inline double CONF_yuef(double x, double y){
	
	return yuef_sept(x, y);
}

/*extern double CONF_yuef(double x, double y);
extern double yuef_sept(double x, double y);*/
extern void CONF_specs(double binx, double biny, double bin2d);
extern void CONF_copycode();
extern char *CONF_prfx(); 

#endif
