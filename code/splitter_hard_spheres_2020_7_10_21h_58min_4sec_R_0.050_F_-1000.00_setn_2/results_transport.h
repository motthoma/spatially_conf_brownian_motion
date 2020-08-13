

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

void calc_transpcoeffs(
		       int setn_per_task, 
		       double t,
		       long int **posshift,
		       long int **negshift,
		       double **posx,
		       double **x_init);

int histogramm_mpi_reduce(int m, 
                          double backshift, 
                          double length, 
                          double bin, 
                          double **positions, 
                          char *fname, 
                          int taskid);

int histogramm2d_mpi_reduce(int m, 
                            double bin2d, 
                            double **positionsx, 
                            double **positionsy,  
                            char *fname, 
                            int taskid);

void print_hist_countercheck(int xcheck, 
			     int ycheck, 
			     int twodcheck, 
			     char *fname_specs);

void TCOEF_copycode(void);
