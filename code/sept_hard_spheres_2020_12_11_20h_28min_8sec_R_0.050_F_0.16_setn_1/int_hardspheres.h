/**ensure that only one header for intra-
 * particle interaction is included in code*/
#ifndef HEADER_INT
#define HEADER_INT

#include "par_sim.h"

#define R_INT (1.0*R_CONF)
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
//extern int intforce(double disti, double dist);
extern void INT_specs(double f_cut);
extern void INT_copycode();
extern char *INT_prfx(); 

#endif
