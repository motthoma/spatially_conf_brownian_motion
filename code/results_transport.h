

/**
 *********************************************************
 *
 * External types and variables
 *
 *********************************************************
 */

/**
 * Structure which contains transport coefficients
 * that are calculated within individual threads.
 **/
typedef struct TAG_TransportCoeffs{
	
        /*average of particle x-positions <x>*/
	long double meanx; 
	/*average of particle speed in x-directions <v>*/
	double meanspeed;
	/*<x^2>*/
	long double meanxsqu; 
	/*<x^3>*/
        long double meanxqub;
        /*mean-squared displacement of x-position (2nd moment): <x^2> - <x>^2*/
        long double msd;
	/*third cumulant of x-position: <x^3> - 3*<x>*<x^2> + 2*<x>^3*/
	long double thirdcum;
        /*effective diffusion coeff.: (<x^2> - <x>^2)/(2*t*B/R)*/
	double deff; 
	/*non-linear mobility mu = <v>/F*/
	double mu;

} T_TransportCoeffs;

extern T_TransportCoeffs tcoeff;

/**
 *********************************************************
 *
 * External functions
 *
 *********************************************************
 */

extern void TCOEF_init();
