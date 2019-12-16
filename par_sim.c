#include <stdio.h>
#include<stdlib.h>
#include <math.h>
/*#include "par_sim.h" */


/* Define parameters for simulation */
int N = 200;
unsigned int setnumb = 1;
int numbtest = 20;
double n = 2*1e5;
int simlong = 20; 
double accur = 1;
double deffaccur = 1e7;
double Fswitch = -3e4;			
double initwidth = 1;
double F = 1;

double time_step(int setnumb, double min_width, double r_core, double force){
double tstepscale, tstep;

if (setnumb == 1) tstepscale = min_width;
else tstepscale = r_core;

if(fabs(force) <= 1/min_width) tstep = tstepscale*tstepscale*pow(10,-3);
else tstep = (tstepscale/fabs(force))*pow(10,-3);

return(tstep);

}

/**struct with parameters used in main function of simulation*/
struct par_specs{
	int N; 
	int setn; 
	double accur; 
	double deffaccur; 
	long long numbtest; 
	double n; 
	double dt; 
	long long testab; 
	long long min_n; 
	double eq_time; 
	double tot_time; 
	int readout_n; 
	double check_time; 
	double readout_time; 
	int numtasks; 
	double f; 
	double initwidth;
};

/**function that provides pointers to parameters used for simulation*/
struct par_specs *par(double n, double dt, int numtasks, int testab, int plotpoints){
	struct par_specs *t_pars;
        t_pars = malloc(sizeof(struct par_specs));
	t_pars->N = N; 
	t_pars->setn = setnumb; 
	t_pars->accur = accur; 
	t_pars->deffaccur = deffaccur; 
	t_pars->numbtest = numbtest; 
	t_pars->n = n; 
	t_pars->dt = dt; 
	t_pars->testab = testab; 
	t_pars->min_n = n+testab*numbtest; 
	t_pars->eq_time = n*dt; 
	t_pars->tot_time = (n+testab*numbtest)*dt;
 	t_pars->readout_n = plotpoints; 
	t_pars->check_time = testab*dt;
	t_pars->readout_time = plotpoints*dt; 
	t_pars->numtasks = numtasks; 
	t_pars->f = F; 
	t_pars->initwidth = initwidth;


return t_pars;
}

/**function used to read simulation parameters in specs file*/
void specs_basic(struct par_specs *t_pars, char *fnamespec){
int a, b, e, h, l, o;
double c, d, f, g, i, j, k, m, n, p, q;

	a = t_pars->N;
        b=t_pars->setn;
        c=t_pars->accur;     
        d=t_pars->deffaccur;
        e=t_pars->numbtest;  
        f=t_pars->n;  
        g=t_pars->dt;  
	h=t_pars->testab;
	i=t_pars->min_n; 
        j=t_pars->eq_time;    
        k=t_pars->tot_time; 
 	l=t_pars->readout_n; 
	m=t_pars->check_time;
	n=t_pars->readout_time;
        o=t_pars->numtasks;
        p=t_pars->f;
        q=t_pars->initwidth;

    FILE *outpspecs;
    outpspecs = fopen(fnamespec, "w");		
    fprintf(outpspecs, "\n\n\nCode can be compiled with:\n\nmpicc -Wall muovertintparallel.c par_sim.c conf_NAME.c int_NAME.c -lgsl -lgslcblas -o muovertintparallel\n\n'masterinteract.py' can be used to create several run-files and start the jobs\n\n\nSimulation Parameters:\n\n");
    

    fprintf(outpspecs, "# of particles: %d\n# of particles per set: %d\n\nAccuracy of mobility: %lf\nAccuracy of Deff: %lf\n# of accuracy checks: %d\n\n# of time steps for eq n=%.2e\ntime step size dt=%.4e\nnumber of steps between two tests: %.1e\nminimum number of simulation steps: %.2e\ntime until checks start: %.3lf\nminimum simulation time: %.3lf\nnumber of steps between two readouts: %.1e\ntime between two accuracy checks: %lf\ntime between readout of two points: %lf\n\n# of parallelized tasks: %d.0\napplied force F: %.2lf\ninitwidth: %.2lf\n\n",a,b,c,d,e,f,g,1.0*h,i,j,k,1.0*l,m,n,o,p,q);
    
fclose(outpspecs);
}

/**function that copies module in working directory*/
void copycode_par(){
char copycode[200];

  sprintf(copycode, "cp ../par_sim* ./");
  system(copycode);

}
