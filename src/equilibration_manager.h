#include <stdbool.h>

#ifndef EQUILIBRATION_MANAGER_H
#define EQUILIBRATION_MANAGER_H


/**
 *********************************************************
 *
 * External types and variables
 *
 *********************************************************
 */

/**
 * structure which contains variables used in the
 * equilibration manager which monitor equilibration
 * of system and calculate quantities that serve as
 * stopping criteria for simulation
 */
typedef struct EquManager{
    int equ_counter;
    int test_start_step;
    double mu_old;
    bool PrintRes;
    bool TestRes;
}T_EquManager;
extern T_EquManager EquManager;

/**
 *********************************************************
 *
 * Internal functions
 *
 *********************************************************
 */

void EquiMan_update_equcounter(double mu_current,
                               double accurarcy);

/**
 *********************************************************
 *
 * External functions
 *
 *********************************************************
 */

void EquiMan_init_params(int stepnumb, int testab);

void EquiMan_init_mu_old(int time_step, double mu);

void EquiManager_set_test_flag(int time_step,
                               int stepnumb,
                               int testab);

#endif /* EQUILIBRATION_MANAGER_H */
