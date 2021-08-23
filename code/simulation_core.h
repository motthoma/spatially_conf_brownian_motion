
#include "comp_gen_header.h"
#include "results_transport.h"
#include "print_routines.h"
/*#include "par_sim.h"
#include "code_handling.h"

#include <stdio.h>
#include <time.h>
#include <stdbool.h>*/
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/**
 *********************************************************
 *
 * External types and variables
 *
 *********************************************************
 */


/**
 *********************************************************
 *
 * Internal functions
 *
 *********************************************************
 */

double reset_pos_time(int setn_per_task, 
                      long int **posshift, 
	    	      long int **negshift);

void adapt_posshifts(int shiftind, 
		     int i, 
		     int j, 
		     long int **posshift, 
		     long int **negshift);

int update_equcounter(double tran_quant, 
		      double tran_quanto, 
		      double accurarcy, 
		      int equcounter);

long int **calloc_2Dlint_array(int m, int n);


/**
 *********************************************************
 *
 * External functions
 *
 *********************************************************
 */

void SIM_init_positions(int setn_per_task, 
                        double **positionx, 
                        double **positiony, 
                        double **xstart, 
                        gsl_rng *r);

void SIM_read_in_positions(int setn_per_task, 
			   double **positionx, 
			   double **positiony, 
			   double **xstart);

void SIM_init_interactions(int setn_per_task, 
                       double **positionx, 
                       double **positiony,
                       double **fintxarray,
                       double **fintyarray);

void SIM_simulation_core(int setn_per_task,
  			 int setn,
			 int taskid, 
                         double **positionx, 
                         double **positiony, 
                         double **xstart, 
                         double **fintxarray,
                         double **fintyarray,
                         gsl_rng *r);

void SIM_copycode(void);

