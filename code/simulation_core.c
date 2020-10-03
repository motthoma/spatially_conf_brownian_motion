#include "par_sim.h"
#include "comp_gen_header.h"
#include "results_transport.h"

#include <stdio.h>
#include <time.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>

void SIM_init_positions(int setn_per_task, 
                       double **positionx, 
                       double **positiony, 
                       double **xstart, 
                       gsl_rng *r)
{
/**
 * function that initializes the particle positions. The particles are uniformly distributed
 * over the whole period length L in a strip of width SimParams.initwidth or the width of
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

  PosValidInit = false;
  for(j = 0; j < setn_per_task; j++){
	  for(kset = 0; kset < SimParams.setnumb; kset++){
		  do{
			  positionx[j][kset] = gsl_rng_uniform(r)*L;
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
	  }
  }	
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
                                          fintxpair = intforce(distx, dist);
                                          fintypair = intforce(disty, dist);
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


