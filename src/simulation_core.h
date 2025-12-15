
#include "comp_gen_header.h"
#include "results_transport.h"
#include "print_routines.h"
/*#include "sim_config.h"
#include "code_handling.h"

#include <stdio.h>
#include <stdbool.h>*/
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#ifndef SIMCOREVALS_H
#define SIMCOREVALS_H

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
typedef struct {
    double **xstart; //x-component of start positions 
} T_EnsembleState;

extern T_EnsembleState EnsembleState;

double reset_pos_time(int setn_per_task, 
                      double **positionx,
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

double **calloc_2Ddouble_array(int m, int n);

/**
 *********************************************************
 *
 * External functions
 *
 *********************************************************
 */

void SIM_alloc_ensemble_state(int setn);

void SIM_init_positions(int setn_per_task, 
                        double **positionx, 
                        double **positiony);

void SIM_read_in_positions(int setn_per_task, 
                           double **positionx, 
                           double **positiony); 

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
                         double **fintxarray,
                         double **fintyarray);

void SIM_copycode(void);

#endif // SIMCOREVALS_H
