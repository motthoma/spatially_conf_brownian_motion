#include "simulation_core.h"
#include "comp_gen_header.h"
/*
#include <stdio.h>
#include <time.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>*/

void SIM_init_positions(int setn_per_task, 
                       double **positionx, 
                       double **positiony, 
                       double **xstart, 
                       gsl_rng *r)
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
  for(j = 0; j < setn_per_task; j++){
	  for(kset = 0; kset < SimParams.setnumb; kset++){
		  do{
			  positionx[j][kset] = gsl_rng_uniform(r)*L_CONF*SimParams.init_max_xpos;
			  positiony[j][kset] = (2*gsl_rng_uniform(r) - 1)*SimParams.initwidth;
			    
			  xo = positionx[j][kset];
			  yo = positiony[j][kset];
		 
			  yue = CONF_yuef(xo, yo);
			 
			  PosValidInit = true; 
			  if(fabs(positiony[j][kset]) >= yue) PosValidInit = false;
			  if((kset > 0) && (PosValidInit == true)){
					    
				  for(ktest = 0; ktest < kset; ktest++){
					  distx = xo - positionx[j][ktest];			
					  disty = yo - positiony[j][ktest];			
					  distinit = sqrt(distx*distx + disty*disty);
					  if(distinit <= 2*R_INT) PosValidInit = false;
				 }
			 }  
			  
		  }while(PosValidInit == false);
		
		  xstart[j][kset] = xo;
		 // printf("x_0:\t%lf\n", xo);
	  }
  }	
  printf("positions initialized!\n"); 
}

void SIM_read_in_positions(int setn_per_task, 
			   double **positionx, 
			   double **positiony, 
			   double **xstart)
{
/**
 * function that reads in positions from file. 
 * File has to be text file with two columns: x-coordinate 
 * is stored in first column, y-coordinate in second column.
 */
  FILE *file;
  int i, j, k;
  double val;
  char val_int;


  file = fopen("test_positions.dat", "r");

  j = 0;
  for(i = 0; i < setn_per_task*SimParams.setnumb; i++)
  {
	for(k = 0; k < 2; k++)
	{
		fscanf(file, "%lf", &val);
		if(k == 0)
		{
		       positionx[i][j] = val;
		}
		else positiony[i][j] = val;
		printf("x: %lf\t y: %lf\n", positionx[i][j], positiony[i][j]);
	}
	j++;
	if(j >= SimParams.setnumb)
	{
		j = 0;
	}
  }

 /* for(i = 0; i < setn_per_task; i++)
  {
	for(j = 0; j < SimParams.setnumb; j++)
	{
		printf("x: %lf\t y: %lf\n", positionx[i][j], positiony[i][j]);
	}
  }*/
  printf("positions read in!\n"); 
}

void SIM_init_interactions(int setn_per_task, 
                       double **positionx, 
                       double **positiony,
                       double **fintxarray,
                       double **fintyarray)
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
                       
  for (j = 0; j < setn_per_task; j++){  
          for(kset = 0; kset < SimParams.setnumb; kset++){
	          fintx = 0;
		  finty = 0;
                  for(ktest = 0; ktest < SimParams.setnumb; ktest++){
                          if(kset != ktest){
                                  distx = positionx[j][kset] - positionx[j][ktest];
                                  disty = positiony[j][kset] - positiony[j][ktest];
                                  dist = sqrt(distx*distx + disty*disty);
                                  if(dist <= INT_CUTOFF){
                                          fintxpair = INT_force(distx, dist);
                                          fintypair = INT_force(disty, dist);
                                          fintx += fintxpair;
                                          finty += fintypair;
                                  }
                          }
                  }
                  fintxarray[j][kset] = fintx;
                  fintyarray[j][kset] = finty;
          }
  }
}

double reset_pos_time(int setn_per_task,
	              double **xstart,	
		      double **xposition,
                      long int **posshift, 
	    	      long int **negshift) 
{
/**
 * Function to reset the simulated system time t and the positions
 * to originate values in order to skip transient effects from
 * start of simulation.
 */
	int i, j;
	double t = 0;
	for(i = 0; i < setn_per_task; i++){
	  for(j = 0; j < SimParams.setnumb; j++){
		 posshift[i][j] = 0; 
		 negshift[i][j] = 0;
		 xstart[i][j] = xposition[i][j];
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

int update_equcounter(double tran_quant, double tran_quanto, double accurarcy, int equcounter)
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

long int **calloc_2Dlint_array(int m, int n){
/**
 * Function that allocates memory for a 2 dimensional array of long ints
 */
long int **array;
int i;

  array = calloc(m, sizeof(double));
  for(i = 0; i < m; i++){
  	array[i] = calloc(n, sizeof(double));
  }

  return array;
}

void SIM_simulation_core(int setn_per_task,
		 	 int setn,
			 int taskid, 
                         double **positionx, 
                         double **positiony, 
                         double **xstart, 
                         double **fintxarray,
                         double **fintyarray,
                         gsl_rng *r)
{
  /* Perform simulation steps until equilibration is reached */

  double u, v, x, y, yo, xo; 
  double xtest, ytest;
  
  double t;
  double dt = SimParams.time_step;  
  long double  sqrt_flucts, f_dt;
  int i;     
  int j;
  int kset; 
  int int_ind; 
  int abb;	
  int abbdeff;	
  double muabbo;	
  double deffabbo;	
  double fintxo;
  double fintyo;
  double fintx;
  double finty;
  double fintxpair;
  double fintypair;
  double yue, distx, disty, dist; 
  int shiftind;
  int ktest;
 
  bool PrintRes;
  bool TestRes;
  bool PosValid;
  
  long int **negshift;
  long int **posshift;
  
  negshift = calloc_2Dlint_array(setn, SimParams.setnumb);
  posshift = calloc_2Dlint_array(setn, SimParams.setnumb);
  
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
	  for (j = 0; j < setn_per_task; j++){	
		  for(kset = 0; kset < SimParams.setnumb; kset++){


			  xo = positionx[j][kset];
			  yo = positiony[j][kset];
			  
			  fintxo = fintxarray[j][kset];
			  fintyo = fintyarray[j][kset];
		  
			 /*
			  *Create random numbers for x- and y-component of noise
			  */  
			  u = gsl_ran_gaussian_ziggurat(r,1);
			  v = gsl_ran_gaussian_ziggurat(r,1);  
			 /* 
			  *Update x- and y-component of position according to
			  *stochastic Euler for Langevin equation
			  */                       
			  x = xo + f_dt + fintxo*dt + sqrt_flucts*u;
			  y = yo + fintyo*dt + sqrt_flucts*v;
			  /*
			   *Calculate value of confinement boundary at current
			   *position x (y-value is needed for non-analytic treatment
			   *of channels with cosine shape
			   */  
			  yue = CONF_yuef(x, y);
			  
			  PosValid = true;

			  /*Check if particle is within effective boundary*/  
			  if (fabs(y) > yue) PosValid = false;
			  /*Check if bottleneck is passed correctly*/	
			  if ((x < 0) && ((fabs(yo) >= BOTTLENECK_WIDTH-R_CONF) || (fabs(y) >= BOTTLENECK_WIDTH-R_CONF))) PosValid = false;
			  if ((x > L_CONF) && ((fabs(yo) >= BOTTLENECK_WIDTH-R_CONF) || (fabs(y) >= BOTTLENECK_WIDTH-R_CONF))) PosValid = false;

			  if(PosValid == true){
				  /*
				   * Adapt position according to period boundary
				   * conditions employed for simulation
				   */
				  shiftind = 0;
				  if(x < 0){
					   shiftind = -1;
					   x += L_CONF;
				  }
				  if(x > L_CONF){
					   shiftind = 1;
					   x -= L_CONF;
				  }
				  
				  /*simulate particle-particle interaction*/ 
				  if (SimParams.setnumb > 1){ 
					    int_ind = 0;
					    fintx = 0;
					    finty = 0;
					  //  for(int_ind = 0; int_ind < SimParams.setnumb; int_ind++){	
					    while((PosValid == true) && (int_ind < SimParams.setnumb)){	
						  xtest = positionx[j][int_ind]; 
						  ytest = positiony[j][int_ind];
 
						  if(int_ind != kset){
							  distx = x - xtest;
							  disty = y - ytest;

							  /* 
							   * search relevant distance according
							   * to minimum image conversion 
							   */
							  if (abs(distx) > 0.5*L_CONF){ 
								distx = distx - L_CONF*(distx/abs(distx));
							  }
							  dist = sqrt(distx*distx + disty*disty);
							  /* 
							   * stay in do-while loop if two particles overlap  
							   * and skip the following break statement
							   */
							  if(dist <= 2*R_INT) PosValid = false;//continue; 
							  
							  if((dist <= INT_CUTOFF) && (PosValid == true)){
								  fintxpair = INT_force(distx, dist);
								  fintypair = INT_force(disty, dist);
								  fintx += fintxpair;
								  finty += fintypair;

							  }
						  }
						  int_ind++;
					    /*close loop over particles of interacting ensemble*/
					    } 
				  }
			  }		   
		          
			  if(shiftind != 0){
				  adapt_posshifts(shiftind, j, kset, posshift, negshift);
                          }
		  
                         /*
                          * Update arrays with positions and inter particle forces
                          */    
		          if(PosValid == true)
			  { 
				  positionx[j][kset] = x; 
				  positiony[j][kset] = y;
				  fintxarray[j][kset] = fintx;
				  fintyarray[j][kset] = finty;
			  }
			  else
			  { 
				  positionx[j][kset] = xo; 
				  positiony[j][kset] = yo;
				  fintxarray[j][kset] = fintxo;
				  fintyarray[j][kset] = fintyo;
			  }

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
		   * */ 
		  RES_calc_transpcoeffs(setn_per_task, 
		   		        t, 
				        posshift,
				        negshift,
				        positionx,
				        xstart);
	  
			 
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
	//	  printf("equcounter updated\n");
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
		  t = reset_pos_time(setn_per_task, 
				     xstart, 
				     positionx, 
				     posshift, 
				     negshift); 
	  } 

 
  /* 
   * Closes while loop over simulation steps if criteria for equilibration are fulfilled
   */
  }while((abb < SimParams.numbtest) || (abbdeff < SimParams.numbtest));
  
  free(negshift);
  free(posshift);
}


void SIM_copycode(){

   char copycode[200];

  sprintf(copycode, "cp ../simulation_core.* ./");
  system(copycode);

}
