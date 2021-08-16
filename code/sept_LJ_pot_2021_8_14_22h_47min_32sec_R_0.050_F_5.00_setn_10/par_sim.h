/*define channel parameters which are needed for all geometries*/
#ifndef HEADER_PAR_SIM
#define HEADER_PAR_SIM

#include <stdbool.h>

/*Parameters that are used for the confinement but also
 * needed for the scaling of the particle-particle
 * interaction and all channel geometries are defined here
 */
/* period length of channel */
#define L_CONF 1.0
/*bottleneck width of periodic channels*/
#define B 0.1*L_CONF

/*radius of particles used in effective spatial boundary */
#define R_CONF (0.5*B)

/*scaling of noise strength depends on ration of bottleneck to particle radius */
#define BOTTRAD ((R_CONF == 0 ? 1.0 : (B/R_CONF)))

/*Order of magnitude for time steps*/
#define TSTEP_BASE 1e-3

/*definition of master thread when mpi parallelization is used*/
#define MASTER 0    

/**
 *********************************************************
 *
 * External types and variables
 *
 *********************************************************
 */
extern unsigned int nbin; 

typedef struct TAG_SimParams{
	int N;
	int setnumb;
	int numbtest;
	int stepnumb;
	int simlong;
	double accur;
	double deffaccur;
	double initwidth;
	double F;
	int plotpoints;
	int testab;
	int reset_stepnumb;
	double tstep_base;
	double time_step;
	int numtasks; 


}T_SimParams;

extern T_SimParams SimParams;

/**
 *********************************************************
 *
 * External functions
 *
 *********************************************************
 */

double PARAMS_time_step(double lscale_conf, double lscale_part);

extern void PARAMS_init(); 
extern bool PARAMS_check_consistency();
extern void PARAMS_basic(char *fnamespec);
extern void PARAMS_copycode();

#endif
