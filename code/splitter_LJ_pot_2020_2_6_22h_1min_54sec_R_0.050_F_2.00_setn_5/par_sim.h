/*define channel parameters which are needed for all geometries*/
#ifndef HEADER_PAR_SIM
#define HEADER_PAR_SIM

//#define L 1.0
#define B 0.1
#define RADF 0.5
#define R_CONF (RADF*B)
#define BOTTRAD ((RADF == 0 ? 1.0 : (1.0/RADF)))

static const double L = 1.0;
extern int N, numbtest, simlong;
extern unsigned int setnumb, nbin; 
extern double n, initwidth, accur, deffaccur, Fswitch, F; 


typedef struct TAG_SimParams{
	int N_test;

}T_SimParams;

extern T_SimParams sim_params;

double time_step(int setnumb, double min_width, double r_core, double force);

struct par_specs{
	int N; 
	int setn; 
	double accur; 
	double deffaccur; 
	int numbtest; 
	double n; 
	double dt; 
	int testab; 
	int min_n; 
	double eq_time; 
	double tot_time; 
	int readout_n; 
	double check_time; 
	double readout_time; 
	int numtasks; 
	double f; 
	double initwidth;
};

struct par_specs *par(double n, double dt, int numtasks, int testab, int plotpoints); 

extern void specs_basic(struct par_specs *t_pars, char *fnamespec);
extern void copycode_par();

#endif
