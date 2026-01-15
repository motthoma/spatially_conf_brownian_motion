#include "comp_gen_header.h"
#include "results_transport.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#ifndef SIMCOREVALS_H
#define SIMCOREVALS_H

/**
 *********************************************************
 *
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
    double **xstart; 
    double **positionx;
    double **positiony;
    double **fintxarray;
    double **fintyarray;
} T_EnsembleState;




void SIM_calculate_inter_particle_forces(const T_SimParams *SimParams, T_EnsembleState *EnsembleState, int j,
                                     int kset,
                                     double x,
                                     double y,
                                     double *fintx,
                                     double *finty);
/**
 *
 *
 *********************************************************
 *
 * External functions
 *
 *********************************************************
 */

T_EnsembleState SIM_alloc_ensemble_state(const T_SimParams *SimParams);

void SIM_init_positions(const T_SimParams *SimParams, T_EnsembleState *EnsembleState); 

void SIM_read_in_positions(const T_SimParams *SimParams, T_EnsembleState *EnsembleState); 

void SIM_init_interactions(const T_SimParams *SimParams, T_EnsembleState *EnsembleState); 

void SIM_simulation_core(const T_SimParams *SimParams, T_EnsembleState *EnsembleState, int taskid); 

void SIM_copycode(void);

#endif 
