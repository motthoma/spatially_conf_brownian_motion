/**ensure that only one header for spatial
 * confinement is included in code*/
#ifndef HEADER_CONF
#define HEADER_CONF

#include <math.h>
#include "sim_config.h"

#define SLOPE 0.9
#define MAX_HALF_WIDTH (SLOPE*L_CONF + BOTTLENECK_WIDTH)

#define SQRT_SHIFT R_CONF*sqrt(1 + SLOPE*SLOPE)
#define R_CONF_SQ  R_CONF*R_CONF
#define Lp (L_CONF - (R_CONF*SLOPE)/(sqrt(1 + SLOPE*SLOPE)))

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

/**function to provide effective boundary of confinement with saw-tooth profile
 * used for the entropic splitter (see Motz et al. J. Chem. Phys. 2014).
 * Effective boundary can only be represented piece wise: In the vicinity of the
 * bottleneck, yueff is a circle with radius R_CONF. In the other sections, we
 * have shifted straight lines parallel to the originial boundary.
 */
static inline double yuef_splitter(double x, double y){

    /*If particle is in central cylinder of width 2B,
     * an evaluation of the eff. boundary is not necessary.
     * Provide value of confinement that ensures y < yueff.
     */
    if(fabs(y) < BOTTLENECK_WIDTH - R_CONF){
        return(MAX_HALF_WIDTH);
    }
    /*evaluate effective boundary*/
    if ((R_CONF < x) && (x < Lp)){
        return (BOTTLENECK_WIDTH + SLOPE*(L_CONF - x) - SQRT_SHIFT);
    }
        
    else if ((0 <= x) && (x <= R_CONF)){
        return (BOTTLENECK_WIDTH - sqrt(R_CONF_SQ - x*x));
    }

    else if ((Lp <= x) && (x <= L_CONF)) {
        return (BOTTLENECK_WIDTH - sqrt(R_CONF_SQ - (x - L_CONF)*(x - L_CONF)));
    }
    else{
        return(5*MAX_HALF_WIDTH);
    }

}

/**
 *********************************************************
 *
 * External functions
 *
 *********************************************************
 */

/*function wrapper that provides generic interface to main*/
static inline double CONF_yuef(double x, double y){
	return yuef_splitter(x, y);
}

/*extern double CONF_yuef(double x, double y);
extern double yuef_splitter(double x, double y);*/
extern void CONF_specs(char *file_confparams);
extern void CONF_copycode();
extern char *CONF_prfx();

#endif
