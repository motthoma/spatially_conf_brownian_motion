#include "simulation_core.h"
#include "comp_gen_header.h"
#include "code_handling.h"
#include "random_numb_gen.h"
#include "sim_config.h"
#include "array_utils.h"
#include <gsl/gsl_rng.h>

T_EnsembleState EnsembleState;

void SIM_alloc_ensemble_state(){
  EnsembleState.xstart = UTILS_calloc_2Ddouble_array(SimParams.n_interact_sets, SimParams.parts_per_set);
  EnsembleState.positionx = UTILS_calloc_2Ddouble_array(SimParams.n_interact_sets, SimParams.parts_per_set);
  EnsembleState.positiony = UTILS_calloc_2Ddouble_array(SimParams.n_interact_sets, SimParams.parts_per_set);
  EnsembleState.fintxarray = UTILS_calloc_2Ddouble_array(SimParams.n_interact_sets, SimParams.parts_per_set);
  EnsembleState.fintyarray = UTILS_calloc_2Ddouble_array(SimParams.n_interact_sets, SimParams.parts_per_set);
}


void SIM_init_positions() 
{
/**
 * function that initializes the particle positions. The particles are uniformly distributed
 * over the whole period length L_CONF in a strip of width SimParams.initwidth or the width of
 * the channel if it is smaller. For interacting particles with hard core radius R_int, 
 * configurations which involve an overlap of two particles are rejected. 
 */
  int j;
  int kset;
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
  for(j = 0; j < SimParams.setn_per_task; j++){
	  for(kset = 0; kset < SimParams.parts_per_set; kset++){
		  do{
               
              EnsembleState.positionx[j][kset] = RNG_get_uniform()*L_CONF*SimParams.init_max_xpos;
			  EnsembleState.positiony[j][kset] = (2*RNG_get_uniform() - 1)*SimParams.initwidth;
			    
			  xo = EnsembleState.positionx[j][kset];
			  yo = EnsembleState.positiony[j][kset];
		 
			  yue = CONF_yuef(xo, yo);
			 
			  PosValidInit = true; 
			  if(fabs(EnsembleState.positiony[j][kset]) >= yue) PosValidInit = false;
			  if((kset > 0) && (PosValidInit == true)){
					    
				  for(ktest = 0; ktest < kset; ktest++){
					  distx = xo - EnsembleState.positionx[j][ktest];			
					  disty = yo - EnsembleState.positiony[j][ktest];			
					  distinit = sqrt(distx*distx + disty*disty);
					  if(distinit <= 2*R_INT) PosValidInit = false;
				 }
			 }  
			  
		  }while(PosValidInit == false);
		
		  EnsembleState.xstart[j][kset] = xo;
		 // printf("x_0:\t%lf\n", xo);
	  }
  }	
  printf("positions initialized!\n"); 
}

void SIM_read_in_positions() 
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
  for(i = 0; i < SimParams.setn_per_task*SimParams.parts_per_set; i++)
  {
	for(k = 0; k < 2; k++)
	{
		fscanf(file, "%lf", &val);
		if(k == 0)
		{
		       EnsembleState.positionx[i][j] = val;
		}
		else EnsembleState.positiony[i][j] = val;
		printf("x: %lf\t y: %lf\n",
               EnsembleState.positionx[i][j],
               EnsembleState.positiony[i][j]);
	}
	j++;
	if(j >= SimParams.parts_per_set)
	{
		j = 0;
	}
  }
  printf("positions read in!\n"); 
}

void SIM_init_interactions() 
{
/**
 * Function that calculates depending on the positions of the particles
 * the initial particle-particle interaction force if such a force e.g.
 * given by a Lennard-Jones interaction is given.
 */
  int j;
  int kset;
  int ktest;
  double distx;
  double disty;
  double dist;
  double fintxpair;
  double fintypair;
  double fintx;
  double finty;
                       
  for (j = 0; j < SimParams.setn_per_task; j++){  
      for(kset = 0; kset < SimParams.parts_per_set; kset++){
          fintx = 0;
          finty = 0;
          for(ktest = 0; ktest < SimParams.parts_per_set; ktest++){
              if(kset != ktest){
                  distx = EnsembleState.positionx[j][kset] - EnsembleState.positionx[j][ktest];
                  disty = EnsembleState.positiony[j][kset] - EnsembleState.positiony[j][ktest];
                  dist = sqrt(distx*distx + disty*disty);
                  if(dist <= INT_CUTOFF){
                      fintxpair = INT_force(distx, dist);
                      fintypair = INT_force(disty, dist);
                      fintx += fintxpair;
                      finty += fintypair;
                  }
              }
          }
          EnsembleState.fintxarray[j][kset] = fintx;
          EnsembleState.fintyarray[j][kset] = finty;
      }
  }
}

double reset_pos_time(long int **posshift, 
                      long int **negshift) 
{
/**
 * Function to reset the simulated system time t and the positions
 * to originate values in order to skip transient effects from
 * start of simulation.
 */
	int i, j;
	double t = 0;
	for(i = 0; i < SimParams.setn_per_task; i++){
	  for(j = 0; j < SimParams.parts_per_set; j++){
		 posshift[i][j] = 0; 
		 negshift[i][j] = 0;
		 EnsembleState.xstart[i][j] = EnsembleState.positionx[i][j];
	  } 
	}
	return t;
}

void adapt_posshifts(int shiftind, 
                     int i, 
                     int j, 
                     long int **posshift, 
                     long int **negshift)
{
/**
 * Function to update shifts which are monitored to calculate 
 * absolute position in x-direction from simulation with cyclic
 * boundary conditions
 */
	/*
	 *shift particle in positive direction if 
	 *position is 'left' of considered channel period   
	 */
	if(shiftind < 0){
	  posshift[i][j]++; 
	} 
	/*
	 *shift particle in negative direction if 
	 *position is 'right' of considered channel period   
	*/
	if(shiftind > 0){
	  negshift[i][j]++;
	}
}

int update_equcounter(double tran_quant,
                      double tran_quanto,
                      double accurarcy,
                      int equcounter)
{
/**
 * Function for the update of the counter that is used to monitore the
 * equilibration of the system. A transport quantity tran_quant is
 * compared with a previous value. If the difference is  below the    
 * demanded accurarcy, the value of the counter is increased.
 * If the difference is larger than twice the accurarcy the counter
 * is decreased
*/ 
	if(fabs(tran_quant - tran_quanto) <= accurarcy){
	  equcounter++;
	}

	if(fabs(tran_quant - tran_quanto) >= accurarcy){
	  equcounter--;
	  if(equcounter < 0) equcounter = 0;
	}
	return equcounter;
}

double perform_sim_step(double old_pos,
                        double force_ext_dt,
                        double force_int_dt,
                        double sqrt_flucts){
    double u, new_pos; 
    /*
     *Create random number for simulation of noise
     */  
     u = RNG_get_gaussian(0.0, 1.0);
     /* 
      *Update new value of position according to
      *stochastic Euler for Langevin equation
      */                       
     new_pos = old_pos + force_ext_dt + force_int_dt + sqrt_flucts*u;

    return new_pos;
}

bool check_pos_validity(double x, double y, double yo){
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

void update_ensemble_state(int j, int kset, double x, double y, double fintx, double finty){ 
    /*
     * Update arrays with positions and inter particle forces
     */     
     EnsembleState.positionx[j][kset] = x; 
     EnsembleState.positiony[j][kset] = y;
     EnsembleState.fintxarray[j][kset] = fintx;
     EnsembleState.fintyarray[j][kset] = finty;
}

bool calculate_inter_particle_forces(int j,
                                     int kset,
                                     double x,
                                     double y,
                                     double *fintx,
                                     double *finty){
    /*function that calculates the intra-particle forces*/
    double distx, disty, dist, xtest, ytest;
    int int_ind = 0;
    bool PosValid = true;
    *fintx = 0;
    *finty = 0;
    while((PosValid == true) && (int_ind < SimParams.parts_per_set)){	
      xtest = EnsembleState.positionx[j][int_ind]; 
      ytest = EnsembleState.positiony[j][int_ind];

      if(int_ind != kset){
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
          /* 
           * stay in do-while loop if two particles overlap  
           * and skip the following break statement
           */
          if(dist <= 2*R_INT){ 
              PosValid = false;
              /* if continue is used, stays in while loop 
               * and search for a valid position is continued*/
              //continue; 
          }
          
          if((dist <= INT_CUTOFF) && (PosValid == true)){
              fintxpair = INT_force(distx, dist);
              fintypair = INT_force(disty, dist);
              *fintx += fintxpair;
              *finty += fintypair;
          }
      }
      int_ind++;
    /*close loop over particles of interacting ensemble*/
    } 

    return PosValid;
}

void shift_pos_for_periodic_bc(double *x, int *shiftind){
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

/* Perform simulation steps until equilibration is reached */
void SIM_simulation_core(int taskid){ 
  double x, y, yo, xo; 
  double t;
  double dt = SimParams.time_step;  
  long double  sqrt_flucts, f_dt;
  int i;     
  int abb;	
  int abbdeff;	
  double muabbo;	
  double deffabbo;	
  double fintx, finty;
  int shiftind;
 
  bool PrintRes;
  bool TestRes;
  bool PosValid;
  
  long int **negshift;
  long int **posshift;
  
  negshift = UTILS_calloc_2Dlint_array(SimParams.n_interact_sets,
                                       SimParams.parts_per_set);
  posshift = UTILS_calloc_2Dlint_array(SimParams.n_interact_sets,
                                       SimParams.parts_per_set);
  
  sqrt_flucts = sqrt(2*BOTTRAD*dt);
  f_dt = SimParams.F*dt;
  
  /*
   * loop over simulation steps.
   * loop is stopped when criterion for equilibration is fulfilled.
   */
  t = 0;
  i = 1;     
  abb = 0;	
  abbdeff = 0;	
  muabbo = 0;	
  deffabbo = 0;	

  printf("start to propagate particles\n");
  do{	
	  t += dt;
	  i++;
	  /* Loop over trajectories */
	  for (int j = 0; j < SimParams.setn_per_task; j++){	
		  for(int kset = 0; kset < SimParams.parts_per_set; kset++){

			  xo = EnsembleState.positionx[j][kset];
			  yo = EnsembleState.positiony[j][kset];
			  
			  fintx = EnsembleState.fintxarray[j][kset];
			  finty = EnsembleState.fintyarray[j][kset];
		  
		 	  /* 
 			   * Perform simulation steps until valid step where no particles overlap
 			   * and particles are within channel is obtained 
 			   */
			  do{
                  x = perform_sim_step(xo, f_dt, fintx*dt, sqrt_flucts);
                  y = perform_sim_step(yo, 0, finty*dt, sqrt_flucts);
                  
                  PosValid = check_pos_validity(x, y, yo);

				  if(PosValid == true){
                      /*
                       * shift positions if x is outside of one period of
                       * channel to account for periodic boundary conditions
                       */ 
                      shift_pos_for_periodic_bc(&x, &shiftind);
					  /* calc interactions here */
					  /* simulate particle-particle interaction */
					  if (SimParams.parts_per_set > 1){ 
                        PosValid = calculate_inter_particle_forces(j, kset, x, y, &fintx, &finty);
					  }
				  }		   
			  }while(PosValid == false);
		          
			  if(shiftind != 0){
                  adapt_posshifts(shiftind, j, kset, posshift, negshift);
              }
              update_ensemble_state(j, kset, x, y, fintx, finty);
		}
 	   /*
        * Close loop over trajectories
        */
	  }
	
      /*
	   *Initialize reference values of mobility and diffusion coefficients
	   *for later judgement of equilibration process.
	   */
	  if((i > SimParams.stepnumb - SimParams.testab) && (i <= SimParams.stepnumb - SimParams.testab + 1)){ 
		  muabbo = tcoeff.mu;
		  deffabbo = tcoeff.deff;
	  }
	
      /*
       * Test progress of equilibration and plot results at certain
       * simulation steps i
       */
      PrintRes = false;
	  if(i % SimParams.plotpoints == 0){
	  	PrintRes = true;
	  }

	  TestRes = false;
	  if((i > SimParams.stepnumb) && (i % SimParams.testab == 0)){
     		TestRes = true;
      }

	  if((TestRes == true) || PrintRes == true){ 

		  /*
		   * call function for calculation of transport coefficients 
		   * such as mobility or mean-squared displacement 
		   */ 
		  RES_calc_transpcoeffs(t, 
                                posshift,
                                negshift,
                                EnsembleState.positionx,
                                EnsembleState.xstart);
          
          /*
           * Update of the equilibration counter that are used to monitore the
           * equilibration of the mobility and diffusivity
           */  
		  if(TestRes == true){ 
 			  abb = update_equcounter(tcoeff.mu, muabbo, SimParams.accur, abb);

			  if(tcoeff.deff > 1.0){  
				  
 			  	  abbdeff = update_equcounter(tcoeff.deff, deffabbo, SimParams.deffaccur*deffabbo, abbdeff);
			  }
			  else{ 
 			  	  abbdeff = update_equcounter(tcoeff.deff, deffabbo, SimParams.deffaccur, abbdeff);
			  }		
			  muabbo = tcoeff.mu;    
			  deffabbo = tcoeff.deff;
          }
 		  /*
           * Plot results to check progress of equilibration 
           */
		  if((PrintRes == true) && (taskid == MASTER)){
			  PRINT_results_over_time(t, 
                                      abb, 
                                      abbdeff);

		  }
	  }

	  /*
       *  Reset position and time information to truncate transient effects from small times 
       */
	  if(i == SimParams.reset_stepnumb){
		  t = reset_pos_time(posshift, 
                             negshift); 
	  } 
  /* 
   * Closes while loop over simulation steps if criteria for equilibration are fulfilled
   */
  }while((abb < SimParams.numbtest) || (abbdeff < SimParams.numbtest));
  
  free(negshift);
  free(posshift);
  free(EnsembleState.xstart);
}

void SIM_copycode(){
    CODEHAND_copy_file_to_dest("simulation_core.c");
    CODEHAND_copy_file_to_dest("simulation_core.h");
}
