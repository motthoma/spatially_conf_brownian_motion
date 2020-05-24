/**ensure that only one header for spatial
 * confinement is included in code*/
#ifndef HEADER_CONF
#define HEADER_CONF

#include <math.h>
#include "par_sim.h"

#define M 0.9
#define MAX_HALF_WIDTH (M*L + B)

#define SQRT_SHIFT R_CONF*sqrt(1 + M*M)
#define R_CONF_SQ  R_CONF*R_CONF
#define Lp (L - (R_CONF*M)/(sqrt(1 + M*M)))

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
if(fabs(y) < B - R_CONF){
	return(MAX_HALF_WIDTH);
}
/*evaluate effective boundary*/
if ((R_CONF <= x) && (x < Lp)){
	return (B + M*(L - x) - SQRT_SHIFT);
}
	
        if (x <= R_CONF){
		return (B - sqrt(R_CONF_SQ - x*x));
	}

		if ((Lp <= x) && (x <= L)) {
			return (B - sqrt(R_CONF_SQ - (x - L)*(x - L)));
		}

return(10);
}

/*function wrapper that provides generic interface to main*/
static inline double yuef_ext(double x, double y){
	
	return yuef_splitter(x, y);
}

/*extern double yuef_ext(double x, double y);
extern double yuef_splitter(double x, double y);*/
extern void specs_conf(double binx, double biny, double bin2d);
extern void copycode_conf();
extern char *prfx_conf();

#endif
