

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
 * Structure which contains parameters 
 * for the histograms used.
 **/
typedef struct TAG_HistogramParams{
	
	/*
	 * width of bins used for spatial discretization for position histograms   
	 */
	
	/* spatial discretization in x-direction (along transport direction) */	
	double binx;
	/* spatial discretization in y-direction (perpendicular transport direction) */	
	double biny;
	/* 2-dimensional spatial discretization */	
	double bin2d; 

} T_HistParams;

extern T_HistParams histparams;

/**
 *********************************************************
 *
 * External functions
 *
 *********************************************************
 */

extern void RES_init();

void RES_calc_transpcoeffs(int setn_per_task, 
		           double t,
		           long int **posshift,
		           long int **negshift,
		           double **posx,
		           double **x_init);

int RES_histogramm_mpi_reduce(int m, 
                              double backshift, 
                              double length, 
                              double bin, 
                              double **positions, 
                              char *fname, 
                              int taskid);

int RES_histogramm2d_mpi_reduce(int m, 
                                double bin2d, 
                                double **positionsx, 
                                double **positionsy,  
                                char *fname, 
                                int taskid);

void RES_print_countercheck(int xcheck, 
			    int ycheck, 
			    int twodcheck, 
			    char *fname_specs);

void RES_copycode(void);
