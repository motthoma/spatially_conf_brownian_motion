/*define channel parameters which are needed for all geometries*/
#ifndef HEADER_PAR_SIM
#define HEADER_PAR_SIM

#include <stdbool.h>

//#define L 1.0
#define B 0.1
#define RADF 0.5
#define R_CONF (RADF*B)
#define BOTTRAD ((RADF == 0 ? 1.0 : (1.0/RADF)))

static const double L = 1.0;
extern unsigned int nbin; 

typedef struct TAG_SimParams{
	int N;
	int setnumb;
	int numbtest;
	int stepnumb;
	int simlong;
	double accur;
	double deffaccur;
	double initwidth;
	double F;
	int plotpoints;
	int testab;
	int reset_stepnumb;
	double time_step;
	int numtasks; 


}T_SimParams;

extern T_SimParams SimParams;

double time_step(double lscale_conf, double lscale_part);

extern void init_simparams(); 
extern bool check_parameter_consistency(); 
extern void specs_basic(char *fnamespec);
extern void copycode_par();

#endif
