#include "simulation_core.h"
#include "comp_gen_header.h"
#include "code_handling.h"
#include "random_numb_gen.h"
#include "sim_config.h"
#include "array_utils.h"
#include "equilibration_manager.h"
#include "print_routines.h"
#include <gsl/gsl_rng.h>

/**
 * Function that allocates memory for the arrays in which the positions and the inter-particle forces are stored.
 * 
 * The arrays are allocated as 2D arrays with dimensions [n_interact_sets][parts_per_set] to allow for a 
 * clear assignment of the particles to the different sets which are simulated in parallel. 
 */
T_EnsembleState SIM_alloc_ensemble_state(const T_SimParams *SimParams){
  T_EnsembleState EnsembleState;
  EnsembleState.xstart = UTILS_calloc_2Ddouble_array(
          SimParams->n_interact_sets,
          SimParams->parts_per_set
          );
  EnsembleState.positionx = UTILS_calloc_2Ddouble_array(
          SimParams->n_interact_sets,
          SimParams->parts_per_set
          );
  EnsembleState.positiony = UTILS_calloc_2Ddouble_array(
          SimParams->n_interact_sets,
          SimParams->parts_per_set
          );
  EnsembleState.fintxarray = UTILS_calloc_2Ddouble_array(
          SimParams->n_interact_sets,
          SimParams->parts_per_set
          );
  EnsembleState.fintyarray = UTILS_calloc_2Ddouble_array(
          SimParams->n_interact_sets,
          SimParams->parts_per_set
          );
  return EnsembleState;
}


/**
 * function that initializes the particle positions. The particles are uniformly distributed
 * over the whole period length L_CONF in a strip of width SimParams->initwidth or the width of
 * the channel if it is smaller. For interacting particles with hard core radius R_int, 
 * configurations which involve an overlap of two particles are rejected. 
 */
void SIM_init_positions(const T_SimParams *SimParams,
                        T_EnsembleState *EnsembleState) 
{
  int set_idx;
  int p_in_set;
  int ktest;
  double xo;
  double yo;
  double yue;
  double distx;
  double disty;
  double distinit;
  bool PosValidInit;

  printf("start to init positions!\n"); 
  
  PosValidInit = false;
  for(set_idx = 0; set_idx < SimParams->setn_per_task; set_idx++){
	  for(p_in_set = 0; p_in_set < SimParams->parts_per_set; p_in_set++){
		  do{
               
              EnsembleState->positionx[set_idx][p_in_set] = RNG_get_uniform()*L_CONF*SimParams->init_max_xpos;
			  EnsembleState->positiony[set_idx][p_in_set] = (2*RNG_get_uniform() - 1)*SimParams->initwidth;
			    
			  xo = EnsembleState->positionx[set_idx][p_in_set];
			  yo = EnsembleState->positiony[set_idx][p_in_set];
		 
			  yue = CONF_yuef(xo, yo);
			 
			  PosValidInit = true; 
			  if(fabs(EnsembleState->positiony[set_idx][p_in_set]) >= yue) PosValidInit = false;
			  if((p_in_set > 0) && (PosValidInit == true)){
					    
				  for(ktest = 0; ktest < p_in_set; ktest++){
					  distx = xo - EnsembleState->positionx[set_idx][ktest];			
					  disty = yo - EnsembleState->positiony[set_idx][ktest];			
					  distinit = sqrt(distx*distx + disty*disty);
					  if(distinit <= 2*R_INT) PosValidInit = false;
				 }
			 }  
			  
		  }while(PosValidInit == false);
		
		  EnsembleState->xstart[set_idx][p_in_set] = xo;
		 // printf("x_0:\t%lf\n", xo);
	  }
  }	
  printf("positions initialized!\n"); 
}

/**
 * function that reads in positions from file. 
 * File has to be text file with two columns: x-coordinate 
 * is stored in first column, y-coordinate in second column.
 */
void SIM_read_in_positions(const T_SimParams *SimParams,
                           T_EnsembleState *EnsembleState) 
{
  FILE *file;
  int i, j, k;
  double val;


  file = fopen("test_positions.dat", "r");

  j = 0;
  for(i = 0; i < SimParams->setn_per_task*SimParams->parts_per_set; i++)
  {
	for(k = 0; k < 2; k++)
	{
		fscanf(file, "%lf", &val);
		if(k == 0)
		{
           EnsembleState->positionx[i][j] = val;
		}
		else EnsembleState->positiony[i][j] = val;
		printf("x: %lf\t y: %lf\n",
               EnsembleState->positionx[i][j],
               EnsembleState->positiony[i][j]);
	}
	j++;
	if(j >= SimParams->parts_per_set)
	{
		j = 0;
	}
  }
  printf("positions read in!\n"); 
}

/**
 * Function to reset the simulated system time t and the positions
 * to originate values in order to skip transient effects from
 * start of simulation.
 */
static double sim_reset_pos_time(const T_SimParams *SimParams,
                                 T_EnsembleState *EnsembleState,
                                 int time_step,
                                 double t,
                                 long int **posshift, 
                                 long int **negshift) 
{
    if(time_step == SimParams->reset_stepnumb){
        int set_idx, p_in_set;
        for(set_idx = 0; set_idx < SimParams->setn_per_task; set_idx++){
          for(p_in_set = 0; p_in_set < SimParams->parts_per_set; p_in_set++){
             posshift[set_idx][p_in_set] = 0; 
             negshift[set_idx][p_in_set] = 0;
             EnsembleState->xstart[set_idx][p_in_set] = EnsembleState->positionx[set_idx][p_in_set];
          } 
        }
        return 0.0;
    }
    else{
        return t;
    }
}

/**
 * Function to update shifts which are monitored to calculate 
 * absolute position in x-direction from simulation with cyclic
 * boundary conditions
 */
static void sim_adapt_posshifts(int shiftind, 
                                int set_idx, 
                                int p_in_set, 
                                long int **posshift, 
                                long int **negshift)
{
	 /*shift particle in positive direction if 
	 *position is 'left' of considered channel period*/
	if(shiftind < 0){
	  posshift[set_idx][p_in_set]++; 
	} 
	
	/*shift particle in negative direction if 
	 *position is 'right' of considered channel period*/
	if(shiftind > 0){
	  negshift[set_idx][p_in_set]++;
	}
}

/**
 * Update arrays with positions and inter particle forces
 */     
static void sim_update_ensemble_state(T_EnsembleState *EnsembleState,
                                      int set_idx,
                                      int p_in_set,
                                      double x,
                                      double y,
                                      double fintx,
                                      double finty){ 
     EnsembleState->positionx[set_idx][p_in_set] = x; 
     EnsembleState->positiony[set_idx][p_in_set] = y;
     EnsembleState->fintxarray[set_idx][p_in_set] = fintx;
     EnsembleState->fintyarray[set_idx][p_in_set] = finty;
}

/**
 * Function that calculates depending on the positions of the particles
 * the initial particle-particle interaction force if such a force e.g.
 * given by a Lennard-Jones interaction is given.
 */
static void sim_calculate_inter_particle_forces(const T_SimParams *SimParams,
                                                T_EnsembleState *EnsembleState) 
{
  int set_idx;
  int p_in_set;
  int ktest;
  double x;
  double y;
  double xtest;
  double ytest;
  double distx;
  double disty;
  double dist;
  double fintxpair;
  double fintypair;
  double fintx;
  double finty;

  for (set_idx = 0; set_idx < SimParams->setn_per_task; set_idx++){  
      for(p_in_set = 0; p_in_set < SimParams->parts_per_set; p_in_set++){
          x =  EnsembleState->positionx[set_idx][p_in_set];
          y =  EnsembleState->positiony[set_idx][p_in_set];
          fintx = 0;
          finty = 0;
          for(ktest = 0; ktest < SimParams->parts_per_set; ktest++){
              xtest =  EnsembleState->positionx[set_idx][ktest];
              ytest =  EnsembleState->positiony[set_idx][ktest];
              if(p_in_set != ktest){
                  distx = x - xtest;
                  disty = y - ytest;

                  /* 
                   * search relevant distance according
                   * to minimum image conversion 
                   */
                  if (fabs(distx) > 0.5*L_CONF){ 
                    distx = distx - L_CONF*(distx/fabs(distx));
                  }
                  dist = sqrt(distx*distx + disty*disty);

                  if(dist <= INT_CUTOFF){
                      fintxpair = INT_force(distx, dist);
                      fintypair = INT_force(disty, dist);
                      fintx += fintxpair;
                      finty += fintypair;
                  }
              }
          }
          EnsembleState->fintxarray[set_idx][p_in_set] = fintx;
          EnsembleState->fintyarray[set_idx][p_in_set] = finty;
      }
  }
}


/**
 * Adapt position according to period boundary
 * conditions employed for simulation
 */
static void sim_shift_pos_for_periodic_bc(double *x, int *shiftind){
    *shiftind = 0;
    if(*x < 0){
       *shiftind = -1;
       *x += L_CONF;
    }
    if(*x > L_CONF){
       *shiftind = 1;
       *x -= L_CONF;
    }
    /* calc interactions here */ 
}

/** 
 * Function to check if a given position of a particle overlaps with any of the other particles in the same set. 
 *
 * This is relevant for the initialization of the positions and for the propagation step if hard-core interactions are present. 
 */
static bool sim_check_particle_overlap(const T_SimParams *SimParams,
                                       T_EnsembleState *EnsembleState,
                                       int set_idx,
                                       int p_in_set_under_check,
                                       double x,
                                       double y) {
    double distx, disty, dist, xtest, ytest;
    for (int p_in_set = 0; p_in_set < SimParams->parts_per_set; p_in_set++) {
        if (p_in_set != p_in_set_under_check) {
            xtest = EnsembleState->positionx[set_idx][p_in_set];
            ytest = EnsembleState->positiony[set_idx][p_in_set];

            distx = x - xtest;
            disty = y - ytest;

            if (fabs(distx) > 0.5 * L_CONF) {
                distx = distx - L_CONF * (distx / fabs(distx));
            }
            dist = sqrt(distx * distx + disty * disty);

            if (dist <= 2 * R_INT) {
                return true; // Overlap detected
            }
        }
    }
    return false; // No overlap
}

/**
 * Function to check if a given position of a particle is valid with respect to the confinement and the bottlenecks. 
 *
 * This is relevant for the propagation step to ensure that particles do not 'tunnel' through the bottlenecks or leave the confinement. 
 */
static bool sim_check_pos_validity(double x, double y, double yo){
    double yue;
    bool PosValid = true;

    /*
     *Calculate value of confinement boundary at current
     *position x (y-value is needed for non-analytic treatment
     *of channels with cosine shape
     */   
    yue = CONF_yuef(x, y);


    /*Check if particle is within effective boundary*/  
    if (fabs(y) > yue) PosValid = false;
    /*Check if bottleneck at x=0 is passed without 'tunneling' through bottleneck*/	
    if ((x < 0) && ((UTILS_max_double(fabs(y), fabs(yo)) >= BOTTLENECK_WIDTH-R_CONF))) PosValid = false;
    /*Check if bottleneck at x=L_CONF is passed without 'tunneling' through bottleneck correctly*/
    if ((x > L_CONF) && ((UTILS_max_double(fabs(y), fabs(yo)) >= BOTTLENECK_WIDTH-R_CONF))) PosValid = false;

    return PosValid;
}

/**
 * Handler function to check if a given position of a particle is valid with respect to the confinement, the bottlenecks and the other particles in the same set.
 */
static bool sim_check_confinement_validity(double x, double y, double yo) {
    return sim_check_pos_validity(x, y, yo);
}

/** 
 * Propagates a particle according to stochastic Euler for Langevin equation with given forces and fluctuations. 
 */
static inline void sim_propagate_particle(double xo,
                                          double f_dt,
                                          double yo,
                                          double fintx_dt,
                                          double finty_dt,
                                          double sqrt_flucts,
                                          double *x,
                                          double *y){
    /*Create random numbers for simulation of 2-dimensional noise*/
    double u_x = RNG_get_gaussian(0.0, 1.0);
    double u_y = RNG_get_gaussian(0.0, 1.0);
    /*Update new value of position according to
     *stochastic Euler for Langevin equation*/
    *x = xo + f_dt + fintx_dt + sqrt_flucts*u_x;
    *y = yo + finty_dt + sqrt_flucts*u_y;
}

/** Propagate a single particle until valid*/
static void sim_perform_valid_step(const T_SimParams *SimParams,
                                   T_EnsembleState *EnsembleState,
                                   int set_idx,
                                   int p_in_set,
                                   double xo,
                                   double yo,
                                   double fintx_dt,
                                   double finty_dt,
                                   double f_dt,
                                   double sqrt_flucts,
                                   double *x_out,
                                   double *y_out,
                                   int *shiftind_out){
    double x, y;
    int shiftind = 0;
    bool PosValid;

    do {
        sim_propagate_particle(xo, f_dt, yo, fintx_dt, finty_dt, sqrt_flucts, &x, &y);
        PosValid = sim_check_confinement_validity(x, y, yo);

        if (PosValid) {
            sim_shift_pos_for_periodic_bc(&x, &shiftind);
            if (SimParams->parts_per_set > 1) {
                if (sim_check_particle_overlap(SimParams, EnsembleState, set_idx, p_in_set, x, y)) {
                    PosValid = false;
                }
            }
        }
    } while (!PosValid);

    *x_out = x;
    *y_out = y;
    *shiftind_out = shiftind;
}

/** Propagate all particles for one set */
static void sim_step_set(const T_SimParams *SimParams,
                         T_EnsembleState *EnsembleState,
                         int set_idx,
                         double f_dt,
                         double sqrt_flucts,
                         long int **posshift,
                         long int **negshift){
    for (int p = 0; p < SimParams->parts_per_set; ++p) {
        double xo = EnsembleState->positionx[set_idx][p];
        double yo = EnsembleState->positiony[set_idx][p];
        double fintx = EnsembleState->fintxarray[set_idx][p];
        double finty = EnsembleState->fintyarray[set_idx][p];
        double x, y;
        int shiftind = 0;

        sim_perform_valid_step(SimParams, EnsembleState, set_idx, p,
                               xo, yo, fintx, finty, f_dt, sqrt_flucts, &x, &y, &shiftind);

        if (shiftind != 0) {
            sim_adapt_posshifts(shiftind, set_idx, p, posshift, negshift);
        }

        sim_update_ensemble_state(EnsembleState, set_idx, p, x, y, fintx, finty);
    }
}

/** 
 * Top-level simulation loop: Contains propagation of all particles 
 * over all time-steps until equilibration 
 */
void SIM_simulation_core(const T_SimParams *SimParams,
                         T_EnsembleState *EnsembleState,
                         int taskid)
{
    long int **posshift = UTILS_calloc_2Dlint_array(SimParams->n_interact_sets, SimParams->parts_per_set);
    long int **negshift = UTILS_calloc_2Dlint_array(SimParams->n_interact_sets, SimParams->parts_per_set);

    long double sqrt_flucts = sqrt(2*BOTTRAD*SimParams->time_step);
    double f_dt = SimParams->F * SimParams->time_step;

    T_EquManager EquManager;
    EQUIMAN_init(&EquManager, SimParams->stepnumb, SimParams->testab);

    double time = 0.0;
    int time_step = 1;

    PRINT_header_for_results_over_time();

    do {
        time += SimParams->time_step;
        ++time_step;

        if (SimParams->parts_per_set > 1) {
            sim_calculate_inter_particle_forces(SimParams, EnsembleState);
        }

        // Propagate each set
        for (int set_idx = 0; set_idx < SimParams->setn_per_task; ++set_idx) {
            sim_step_set(SimParams, EnsembleState, set_idx, f_dt, sqrt_flucts, posshift, negshift);
        }

        EQUIMAN_update_mu_old(&EquManager, time_step, tcoeff.mu);
        PRINT_set_print_flag(time_step);
        EQUIMAN_set_test_flag(&EquManager, time_step, SimParams->stepnumb, SimParams->testab);

        if (EquManager.TestRes || Print.PrintRes) {
            RES_calc_transpcoeffs(time, posshift, negshift, EnsembleState->positionx, EnsembleState->xstart);
        }

        EQUIMAN_update_counter(&EquManager, tcoeff.mu, SimParams->accur);

        if (taskid == MASTER) {
            PRINT_results_over_time(time, EquManager.equ_counter);
        }

        time = sim_reset_pos_time(SimParams, EnsembleState, time_step, time, posshift, negshift);

    } while (EquManager.equ_counter < SimParams->patience);

    UTILS_free_2Dlint_array(negshift);
    UTILS_free_2Dlint_array(posshift);
    UTILS_free_2Ddouble_array(EnsembleState->xstart);
}

/**
 * copies the code specific for the simulation_core to 
 * the destination folder for documentation purposes
 */
void SIM_copycode(){
    CODEHAND_copy_file_to_dest("simulation_core.c");
    CODEHAND_copy_file_to_dest("simulation_core.h");
}
