/**ensure that only one header for intra-
 * particle interaction is included in code*/
#ifndef HEADER_INT
#define HEADER_INT

#include <math.h>
#define R_INT (R_CONF)
#define INT_CUTOFF (0.2*L)
#define EPS_L (1.0)
#define LJMIN (0.5*B)
#define LJMINPOW (pow(LJMIN,6))
#define LJPREFAC 12.0*EPS_L*LJMINPOW



extern double fintx, finty, fintxpair, fintypair;


/*extern double intforce_lj(double disti, double dist);
extern double intforce(double disti, double dist);*/

/**function for the calculation of the intra-particle force arising
 * from a Lennard-Jones Potential. The potential is described by two
 * parameters, position of minimum LJ_MIN and depth of minimum  LJ_EPS.
 * The potential is truncated at INT_CUTOFF. 
 */
static inline double intforce_lj(double dist1d, double dist2d){

	return(LJPREFAC*(LJMINPOW*pow(dist2d,-14)-pow(dist2d,-8))*dist1d);

}

static inline double intforce(double dist1d, double dist2d){

	return intforce_lj(dist1d, dist2d);

}
extern void specs_int(double f_cut);
extern void copycode_int();
extern char *prfx_int(); 

#endif
