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
#define BOTTLENECK_WIDTH 0.1*L_CONF

/*radius of particles used in effective spatial boundary */
#define R_CONF (0.5*BOTTLENECK_WIDTH)

/*scaling of noise strength depends on ration of bottleneck to particle radius */
#define BOTTRAD ((R_CONF == 0 ? 1.0 : (BOTTLENECK_WIDTH/R_CONF)))

/*Order of magnitude for time steps*/
#define TSTEP_BASE 1e-3

/*definition of master thread when mpi parallelization is used*/
#define MASTER 0    

/*path to dir location where dirs with simulation results are stored */
#define RUNS_DIR "../runs"

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
	int parts_per_set;
	int numbtest;
	int stepnumb;
	int simlong;
	double accur;
	double deffaccur;
	double initwidth;
	double init_max_xpos;
	double F;
	int plotpoints;
	int testab;
	int reset_stepnumb;
	double time_step;
	int numtasks; 
    int setn_per_task;
    int n_interact_sets;
}T_SimParams;

extern T_SimParams SimParams;

/**
 *********************************************************
 *
 * External functions
 *
 *********************************************************
 */

double SIMCONFIG_time_step(double lscale_conf, double lscale_part);

extern void SIMCONFIG_init(int tasks); 
extern bool SIMCONFIG_check_consistency();
extern void SIMCONFIG_write_specs();
extern void SIMCONFIG_copy_code();

#endif
