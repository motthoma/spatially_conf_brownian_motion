/**ensure that only one header for intra-
 * particle interaction is included in code*/
#ifndef HEADER_INT
#define HEADER_INT

#include <math.h>
#include "par_sim.h"

#define R_INT (R_CONF)
#define INT_CUTOFF (0.5*L_CONF)
#define EPS_L (0.3)
#define LJMIN (0.5*B)

/* defines that are precomputed for better
 * computational performance */
#define LJMINPOW (pow(LJMIN,6))
#define LJPREFAC 12.0*EPS_L*LJMINPOW
#define FSHIFT LJPREFAC*(LJMINPOW*pow(INT_CUTOFF,-14) - pow(INT_CUTOFF,-8))*INT_CUTOFF

/**
 *********************************************************
 *
 * External types and variables
 *
 *********************************************************
 */

extern double fintx, finty, fintxpair, fintypair;

/**
 *********************************************************
 *
 * Internal functions
 *
 *********************************************************
 */

static inline double intforce_lj(double dist1d, double dist2d){
/** Function that returns value of truncated LJ force of two particles which
 *  that move in two dimensions and are separated by the distances dist1d 
 *  and dist2d in each dimension */

	return(LJPREFAC*(LJMINPOW*pow(dist2d,-14) - pow(dist2d,-8))*dist1d - FSHIFT);

}

/**
 *********************************************************
 *
 * External functions
 *
 *********************************************************
 */


/*extern double intforce_lj(double disti, double dist);
extern double intforce(double disti, double dist);*/

/**function for the calculation of the intra-particle force arising
 * from a Lennard-Jones Potential. The potential is described by two
 * parameters, position of minimum LJ_MIN and depth of minimum  LJ_EPS.
 * The potential is truncated at INT_CUTOFF. 
 */

static inline double INT_force(double dist1d, double dist2d){

	return intforce_lj(dist1d, dist2d);

}

extern void INT_specs(char *intspecs);
extern void INT_copycode_int();
extern char *INT_prfx(); 

#endif
