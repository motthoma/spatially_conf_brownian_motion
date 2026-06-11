#include "simulation_core.h"
#include "comp_gen_header.h"
#include "code_handling.h"
#include "random_numb_gen.h"
#include "sim_config.h"
#include "array_utils.h"
#include "equilibration_manager.h"
#include "print_routines.h"
#include <gsl/gsl_rng.h>

T_EnsembleState SIM_alloc_ensemble_state(const T_SimParams *SimParams){
  T_EnsembleState EnsembleState;
  EnsembleState.xstart = UTILS_calloc_2Ddouble_array(SimParams->n_interact_sets, SimParams->parts_per_set);
  EnsembleState.positionx = UTILS_calloc_2Ddouble_array(SimParams->n_interact_sets, SimParams->parts_per_set);
  EnsembleState.positiony = UTILS_calloc_2Ddouble_array(SimParams->n_interact_sets, SimParams->parts_per_set);
  EnsembleState.fintxarray = UTILS_calloc_2Ddouble_array(SimParams->n_interact_sets, SimParams->parts_per_set);
  EnsembleState.fintyarray = UTILS_calloc_2Ddouble_array(SimParams->n_interact_sets, SimParams->parts_per_set);
  return EnsembleState;
}


void SIM_init_positions(const T_SimParams *SimParams,
                        T_EnsembleState *EnsembleState) 
{
/**
 * function that initializes the particle positions. The particles are uniformly distributed
 * over the whole period length L_CONF in a strip of width SimParams->initwidth or the width of
 * the channel if it is smaller. For interacting particles with hard core radius R_int, 
 * configurations which involve an overlap of two particles are rejected. 
 */
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

void SIM_read_in_positions(const T_SimParams *SimParams,
                           T_EnsembleState *EnsembleState) 
{
/**
 * function that reads in positions from file. 
 * File has to be text file with two columns: x-coordinate 
 * is stored in first column, y-coordinate in second column.
 */
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

void SIM_init_interactions(const T_SimParams *SimParams,
                           T_EnsembleState *EnsembleState) 
{
/**
 * Function that calculates depending on the positions of the particles
 * the initial particle-particle interaction force if such a force e.g.
 * given by a Lennard-Jones interaction is given.
 */
  int set_idx;
  int p_in_set;
  int ktest;
  double distx;
  double disty;
  double dist;
  double fintxpair;
  double fintypair;
  double fintx;
  double finty;
                       
  for (set_idx = 0; set_idx < SimParams->setn_per_task; set_idx++){  
      for(p_in_set = 0; p_in_set < SimParams->parts_per_set; p_in_set++){
          fintx = 0;
          finty = 0;
          for(ktest = 0; ktest < SimParams->parts_per_set; ktest++){
              if(p_in_set != ktest){
                  distx = EnsembleState->positionx[set_idx][p_in_set] - EnsembleState->positionx[set_idx][ktest];
                  disty = EnsembleState->positiony[set_idx][p_in_set] - EnsembleState->positiony[set_idx][ktest];
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

double sim_reset_pos_time(const T_SimParams *SimParams,
                          T_EnsembleState *EnsembleState,
                          int time_step,
                          double t,
                          long int **posshift, 
                          long int **negshift) 
{
/**
 * Function to reset the simulated system time t and the positions
 * to originate values in order to skip transient effects from
 * start of simulation.
 */

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

void sim_adapt_posshifts(int shiftind, 
                         int set_idx, 
                         int p_in_set, 
                         long int **posshift, 
                         long int **negshift)
{
/**
 * Function to update shifts which are monitored to calculate 
 * absolute position in x-direction from simulation with cyclic
 * boundary conditions
 */
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

double sim_perform_sim_step(double old_pos,
                            double force_ext_dt,
                            double force_int_dt,
                            double sqrt_flucts){
    double u, new_pos; 
    /*Create random number for simulation of noise*/  
     u = RNG_get_gaussian(0.0, 1.0);
     /*Update new value of position according to
      *stochastic Euler for Langevin equation*/
     new_pos = old_pos + force_ext_dt + force_int_dt + sqrt_flucts*u;

    return new_pos;
}

bool sim_check_pos_validity(double x, double y, double yo){
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

void sim_update_ensemble_state(T_EnsembleState *EnsembleState,
                               int set_idx,
                               int p_in_set,
                               double x,
                               double y,
                               double fintx,
                               double finty){ 
    /*
     * Update arrays with positions and inter particle forces
     */     
     EnsembleState->positionx[set_idx][p_in_set] = x; 
     EnsembleState->positiony[set_idx][p_in_set] = y;
     EnsembleState->fintxarray[set_idx][p_in_set] = fintx;
     EnsembleState->fintyarray[set_idx][p_in_set] = finty;
}

void SIM_calculate_inter_particle_forces(const T_SimParams *SimParams,
                                         T_EnsembleState *EnsembleState,
                                         int set_idx,
                                         int p_in_set,
                                         double x,
                                         double y,
                                         double *fintx,
                                         double *finty){
    /*function that calculates the intra-particle forces*/
    double distx, disty, dist, xtest, ytest;
    double fintxpair, fintypair;
    *fintx = 0;
    *finty = 0;
    for(int int_ind = 0; int_ind < SimParams->parts_per_set; int_ind++){ 
      xtest = EnsembleState->positionx[set_idx][int_ind]; 
      ytest = EnsembleState->positiony[set_idx][int_ind];

      if(int_ind != p_in_set){
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
              *fintx += fintxpair;
              *finty += fintypair;
          }
      }
    } 
}
void sim_shift_pos_for_periodic_bc(double *x, int *shiftind){
    /*
    * Adapt position according to period boundary
    * conditions employed for simulation
    */
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

static bool sim_check_confinement_validity(double x, double y, double yo) {
    return sim_check_pos_validity(x, y, yo);
}

static void sim_propagate_particle(double xo,
                                   double f_dt,
                                   double yo,
                                   double fintx_dt,
                                   double finty_dt,
                                   double sqrt_flucts,
                                   double *x,
                                   double *y) {
    *x = sim_perform_sim_step(xo, f_dt, fintx_dt, sqrt_flucts);
    *y = sim_perform_sim_step(yo, 0, finty_dt, sqrt_flucts);
}

/* Perform simulation steps until equilibration is reached */
void SIM_simulation_core(const T_SimParams *SimParams,
                         T_EnsembleState *EnsembleState,
                         int taskid){ 
  double x, y, yo, xo; 
  double time;
  double dt = SimParams->time_step;  
  long double  sqrt_flucts, f_dt;
  int time_step;     
  int test_start_step;
  double fintx, finty;
  int shiftind;
 
  T_EquManager EquManager;
 
  bool PosValid;
  
  long int **negshift;
  long int **posshift;
  
  negshift = UTILS_calloc_2Dlint_array(SimParams->n_interact_sets,
                                       SimParams->parts_per_set);
  posshift = UTILS_calloc_2Dlint_array(SimParams->n_interact_sets,
                                       SimParams->parts_per_set);
  
  sqrt_flucts = sqrt(2*BOTTRAD*dt);
  f_dt = SimParams->F*dt;
  
  /*
   * loop over simulation steps.
   * loop is stopped when criterion for equilibration is fulfilled.
   */
  time = 0;
  time_step = 1;
  PRINT_header_for_results_over_time(); 
  EQUIMAN_init(&EquManager, SimParams->stepnumb, SimParams->testab);

  printf("loop over trajectory started.");
  do{
	  time += dt;
	  time_step++;
	  /* Loop over trajectories */
	  for (int set_idx = 0; set_idx < SimParams->setn_per_task; set_idx++){
		  for(int p_in_set = 0; p_in_set < SimParams->parts_per_set; p_in_set++){

			  xo = EnsembleState->positionx[set_idx][p_in_set];
			  yo = EnsembleState->positiony[set_idx][p_in_set];
			  
			  fintx = EnsembleState->fintxarray[set_idx][p_in_set];
			  finty = EnsembleState->fintyarray[set_idx][p_in_set];
		  
			  /* 
 			   * Perform simulation steps until valid step where no particles overlap
 			   * and particles are within channel is obtained 
 			   */
			  do{
                  sim_propagate_particle(xo, f_dt, yo, fintx * dt, finty * dt, sqrt_flucts, &x, &y);

                  PosValid = sim_check_confinement_validity(x, y, yo);

				  if(PosValid == true){
                      sim_shift_pos_for_periodic_bc(&x, &shiftind);
					  if (SimParams->parts_per_set > 1){
                          if(sim_check_particle_overlap(SimParams,
                                                        EnsembleState,
                                                        set_idx,
                                                        p_in_set,
                                                        x,
                                                        y)){
                            PosValid = false;
                          }
					  }
				  }
			  }while(PosValid == false);

			  if (SimParams->parts_per_set > 1){
				SIM_calculate_inter_particle_forces(SimParams,
                                                    EnsembleState,
                                                    set_idx,
                                                    p_in_set,
                                                    x,
                                                    y,
                                                    &fintx,
                                                    &finty);
			  }

			  if(shiftind != 0){
                  sim_adapt_posshifts(shiftind,
                                      set_idx,
                                      p_in_set,
                                      posshift,
                                      negshift);
              }
              sim_update_ensemble_state(EnsembleState,
                                        set_idx,
                                        p_in_set,
                                        x,
                                        y,
                                        fintx,
                                        finty);
		}
 	   /* Close loop over trajectories*/
	  }
      EQUIMAN_update_mu_old(&EquManager, time_step, tcoeff.mu);
 
      PRINT_set_print_flag(time_step);
      EQUIMAN_set_test_flag(&EquManager,
                            time_step,
                            SimParams->stepnumb,
                            SimParams->testab);

	  if((EquManager.TestRes == true) || (Print.PrintRes == true)){ 
		  /* call function for calculation of transport coefficients 
		   * such as mobility or mean-squared displacement */ 
		  RES_calc_transpcoeffs(time, 
                                posshift,
                                negshift,
                                EnsembleState->positionx,
                                EnsembleState->xstart);
	  }

      /* Update of the equilibration counter that are used to monitore the
       * equilibration of the mobility and diffusivity*/
       EQUIMAN_update_counter(&EquManager, tcoeff.mu, SimParams->accur);

      /* Plot results to check progress of equilibration*/
      if(taskid == MASTER){
          PRINT_results_over_time(time, EquManager.equ_counter); 
      }

	  /*  Reset position and time information to truncate
       *  transient effects from small times */
      time = sim_reset_pos_time(SimParams,
                                EnsembleState,
                                time_step,
                                time,
                                posshift, 
                                negshift); 
  /* Closes while loop over simulation steps if criteria for equilibration are fulfilled */
  }while(EquManager.equ_counter < SimParams->patience);
  
  free(negshift);
  free(posshift);
  free(EnsembleState->xstart);
}

void SIM_copycode(){
    CODEHAND_copy_file_to_dest("simulation_core.c");
    CODEHAND_copy_file_to_dest("simulation_core.h");
}
