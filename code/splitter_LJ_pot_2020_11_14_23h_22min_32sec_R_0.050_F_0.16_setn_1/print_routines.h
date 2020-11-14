
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
	/*state index to monitor if header line
         *or ongoing value has to be printed
         */
	int state;
        /*array to store name of .dat file containing mobility and other coeffs*/
	char fname [60];
        /*array to store name of .dat file containing moments of positions*/
	char fnamemom [60];

} printres;


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

void PRINT_results_over_time(double t, 
                             int abb, 
                             int abbdeff);

void PRINT_copycode(void);

