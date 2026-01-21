#ifndef PRINT_ROUTINES_H
#define PRINT_ROUTINES_H

#include <time.h>
#include <stdbool.h>

/**
 *********************************************************
 *
 * External types and variables
 *
 *********************************************************
 */

/**
 * structure which contains variables used in the
 * printing functions which print the resulting
 * transport coefficients in .dat files
 */
struct PrintResults{
    /*array to store name of .dat file containing mobility and other coeffs*/
	char fname [256];
    /*array to store name of .dat file containing moments of positions*/
	char fnamemom [256];
    /*flag that controls if snapshot of transport parameters is printed*/
    bool PrintRes;
};

extern struct PrintResults Print;

/**
 *********************************************************
 *
 * Internal functions
 *
 *********************************************************
 */


/**
 *********************************************************
 *
 * External functions
 *
 *********************************************************
 */


void PRINT_positions(double **posx, double **posy);

void PRINT_muoverf(double muall, double deffall);


void PRINT_set_print_flag(int time_step);

void PRINT_header_for_results_over_time(); 

void PRINT_results_over_time(double t, 
                             int abb); 

void PRINT_runtime(clock_t start);

void PRINT_runtime_threads(clock_t start,
                           int numtasks, 
                           int taskid); 

void PRINT_resallthreads(long double msdall, 
                         double meanspeedall, 
                         double muall, 
                         double deffall, 
                         long double meanxall, 
                         long double meanxsquall, 
                         long double thirdcumall);

void PRINT_copycode(void);

#endif /* PRINT_ROUTINES_H */
