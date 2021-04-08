#include <stdio.h>
#include<stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "par_sim.h" 


T_SimParams SimParams;

/**Function that initializes hard coded parameters of simulation*/
void PARAMS_init(){
	SimParams.N = 100;
	SimParams.numbtest = 20;
	SimParams.stepnumb = 0.5*1e6;
	SimParams.simlong = 20;
	SimParams.accur = 1;
	SimParams.deffaccur = 1e7;
	SimParams.initwidth = 1.0;
	
	SimParams.plotpoints =  SimParams.stepnumb/50; 
	SimParams.testab = SimParams.plotpoints*1;
        SimParams.reset_stepnumb = 10*SimParams.testab;
	
	if(fabs(SimParams.F) <= 0.1) SimParams.stepnumb = SimParams.simlong*SimParams.stepnumb;

	if(SimParams.F <= -2*pow(10,4)) SimParams.stepnumb = 0.5*SimParams.simlong*SimParams.stepnumb;

}

bool PARAMS_check_consistency(){
/**
 * Checks consistency of some parameters
 */
	bool ParaFlag = true;
	if(SimParams.testab % SimParams.plotpoints != 0){ 
	  printf("testab modulo SimParams.plotpoints \n");
	  ParaFlag = false;
	}
	if(SimParams.N % (SimParams.numtasks*SimParams.setnumb) != 0){ 
	  printf("N modulo numtasks*setnumb not zero!\n");
	  ParaFlag = false;
	}
	if(ParaFlag == false){
		printf("Chosen simulation parameters are inconsistent!\n");
		return -1;
	}	

}


/**Function that initializes value of simulation time steps.
 * Time steps are calculated by means of confinement and particle quantitites.
 * Therefore, the function has to be called in main, where header files for 
 * confinement and particle-particle interactions are included. */
double PARAMS_time_step(double lscale_conf, double lscale_part){
double tstepscale, tstep;

	if (SimParams.setnumb == 1) tstepscale = lscale_conf;
	else tstepscale = lscale_part;

	if(fabs(SimParams.F) <= 1/lscale_conf) tstep = tstepscale*tstepscale*pow(10,-3);
	else tstep = (tstepscale/fabs(SimParams.F))*pow(10,-3);

return(tstep*L);
}

/**function used to read simulation parameters in specs file*/
void PARAMS_basic(char *fnamespec){
/*some quantities that are only calculated for informational purpose in specs file,
 * but which are not needed for simulation
 */	
int min_n;
double eq_time;
double tot_time;
double check_time;
double readout_time;

    min_n = SimParams.stepnumb + SimParams.testab*SimParams.numbtest; 
    eq_time = (SimParams.stepnumb - SimParams.reset_stepnumb)*SimParams.time_step; 
    tot_time = (SimParams.stepnumb - SimParams.reset_stepnumb + SimParams.testab*SimParams.numbtest)*SimParams.time_step;
    check_time = SimParams.testab*SimParams.time_step;
    readout_time = SimParams.plotpoints*SimParams.time_step; 

    FILE *outpspecs;
    outpspecs = fopen(fnamespec, "w");		
    fprintf(outpspecs, "\n\n\nCode can be compiled with:\n\nmpicc -Wall muovertintparallel.c par_sim.c conf_NAME.c int_NAME.c -lgsl -lgslcblas -o muovertintparallel\n\n'masterinteract.py' can be used to create several run-files and start the jobs\n\n\nSimulation Parameters:\n\n");


    fprintf(outpspecs, "# of particles: %d\n"
                       "# of particles per set: %d\n\n"
		       "Accuracy of mobility: %lf\n"
		       "Accuracy of Deff: %lf\n\n",
		       SimParams.N,
		       SimParams.setnumb,
		       SimParams.accur,
		       SimParams.deffaccur);


    fprintf(outpspecs, 
       		       "number of steps until equilibration checks start: %.2e\n"
		       "time step size for propagation of Langevin equation: %.2e\n"
		       "time until equilibration is reset: %lf\n"
		       "number of checks to validate equilibration: %d\n"
		       "number of steps between two tests: %.2e\n"
		       "number of steps between two readouts: %.2e\n"
		       "time between two accuracy checks: %lf\n"
		       "time between readout of two points: %lf\n"
		       "minimum number of simulation steps: %.2e\n"
		       "minimum simulation time: %.3lf\n\n",
		       (double) SimParams.stepnumb,
		       SimParams.time_step,
		       eq_time,
		       SimParams.numbtest,
		       (double) SimParams.testab,
		       (double) SimParams.plotpoints,
		       check_time,
		       readout_time,
		       (double) min_n,
		       tot_time);

    fprintf(outpspecs, 
                      "# of parallelized tasks: %d\n"
		       "applied force F: %.2lf\n"
		       "width of strip of initial particle distributions: %.2lf\n\n", 
		       SimParams.numtasks,
		       SimParams.F,
		       SimParams.initwidth);

    
fclose(outpspecs);
}

void PARAMS_copycode(){
/**
 * copies module in working directory
 */
char copycode[200];

  sprintf(copycode, "cp ../par_sim* ./");
  system(copycode);

}
