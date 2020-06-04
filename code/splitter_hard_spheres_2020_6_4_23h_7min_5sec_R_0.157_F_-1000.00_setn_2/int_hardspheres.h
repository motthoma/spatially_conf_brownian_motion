/**ensure that only one header for intra-
 * particle interaction is included in code*/
#ifndef HEADER_INT
#define HEADER_INT

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

static inline double intforce(double disti, double dist){

	return(0);

}
//extern int intforce(double disti, double dist);
extern void specs_int(double f_cut);
extern void copycode_int();
extern char *prfx_int(); 

#endif
