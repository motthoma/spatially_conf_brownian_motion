/*define channel parameters which are needed for all geometries*/
#ifndef HEADER_PAR_SIM
#define HEADER_PAR_SIM

//#define L 1.0
#define B 0.1
#define RADF 0.5
#define R_CONF (RADF*B)
#define BOTTRAD ((RADF == 0 ? 1.0 : (1.0/RADF)))

static const double L = 1.0;
extern unsigned int nbin; 
//extern double initwidth, F; 


typedef struct TAG_SimParams{
	int N;
	int setnumb;
	int numbtest;
	int stepnumb;
	int simlong;
	double accur;
	double deffaccur;
//	double Fswitch;
	double initwidth;
	double F;
	int plotpoints;
	int testab;
	int eq_stepnumb;
	double eq_time;
	double time_step;
	double tot_time;
//	int readout_n;
	double check_time;
	double readout_time; 
	int numtasks; 


}T_SimParams;

extern T_SimParams SimParams;

double time_step(double lscale_conf, double lscale_part);

/*struct par_specs{
//	int N; 
//	int setn; 
//	double accur; 
//	double deffaccur; 
//	int numbtest; 
//	double n; 
//	double dt; 
//	int testab; 
//	int min_n; 
//	double eq_time; 
//	double tot_time; 
//	int readout_n; 
//	double check_time; 
//	double readout_time; 
//	int numtasks; 
//	double f; 
//	double initwidth;
};*/

//struct par_specs *par(); 

extern void init_simparams(); 
extern void specs_basic(char *fnamespec);
extern void copycode_par();

#endif
