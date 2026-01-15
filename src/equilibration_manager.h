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
typedef struct equilibration_manager{
    int equ_counter;
    int test_start_step;
    double mu_old;
    bool PrintRes;
    bool TestRes;
}T_EquManager;


/**
 *********************************************************
 *
 * Internal functions
 *
 *********************************************************
 */

void EQUIMAN_update_counter(T_EquManager *EquManager, double mu_current,
                               double accurarcy);

/**
 *********************************************************
 *
 * External functions
 *
 *********************************************************
 */

void EQUIMAN_init(T_EquManager *EquManager, int stepnumb, int testab);

void EQUIMAN_update_mu_old(T_EquManager *EquManager, int time_step, double mu);

void EQUIMAN_set_test_flag(T_EquManager *EquManager, int time_step,
                               int stepnumb,
                               int testab);

#endif /* EQUILIBRATION_MANAGER_H */
