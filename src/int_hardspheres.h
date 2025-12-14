/**ensure that only one header for intra-
 * particle interaction is included in code*/
#ifndef HEADER_INT
#define HEADER_INT

#include "sim_config.h"

/* particle radius relevant for inter-particle
 * interaction. Usualy the same as the particle
 * radius relevant for the confinement
 */
#define R_INT (1.0*R_CONF)

/* spatial cut-off of interaction to avoid 
 * long-range interactions over infinite 
 * distances. Variable not relavant for 
 * hard spheres.
 */
#define INT_CUTOFF 0

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
 * External functions
 *
 *********************************************************
 */

static inline double INT_force(double disti, double dist){
	return(0);
}

extern void INT_specs(char *intspecs);
extern void INT_copycode();
extern char *INT_prfx(); 

#endif
