/*Calculation of the averagespeed for a brownian partical with finite radius in 2d in a confinement, 
trajektories are calculatet parallel*/


#include "par_sim.h"
//#include "comp_gen_header.h"
//#include "results_transport.h"
#include "code_handling.h"
#include "simulation_core.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
/*#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>*/


void print_runtime(clock_t start)
{
/**
 * prints run time of code into file
 */
 
  int timediff;
  int timediff_all;
  clock_t end = clock();
  
  timediff = (int)((end-start) / CLOCKS_PER_SEC);
  timediff_all = (int)(SimParams.numtasks*(end-start) / CLOCKS_PER_SEC);

  FILE *outpspecs;
  outpspecs=fopen("muovert_specs.dat", "a");
  fprintf(outpspecs, "\nComputing time: %d days %dhours %dmin %dsec\n", timediff/(3600*24), (timediff/3600)%24, 
                                                                        (timediff/60)%60, timediff%60);

  fprintf(outpspecs, "Total computing time of all threads: %d days %dhours %dmin %dsec\n\n", timediff_all/(3600*24), 
                                                                                             (timediff_all/3600)%24, 
                                                                                             (timediff_all/60)%60, 
                                                                                              timediff_all%60);
  fclose(outpspecs);

}

void print_positions(int m, double **posx, double **posy){
/**
 * prints particle positions to file
 */
  int i;
  int j;
  FILE *outpos;
  outpos=fopen("positions.dat", "a");
  fprintf(outpos, "#xpositions\t ypositions\n");

  for(i = 0; i < m; i++){
	for(j = 0; j < SimParams.setnumb; j++){
		fprintf(outpos, "%.5lf\t %.5lf\n", posx[i][j], posy[i][j]);
		
	}

  }
  
  fclose(outpos);

}

void print_runtime_threads(clock_t start,
			   int numtasks, 
			   int taskid) 

{
/**
 * prints individual runtime of each thread
 */
clock_t end;
int timediff;

  if(numtasks > 1){
        /**
         * measure program run time of task
         */ 
	end = clock();
	timediff = (int)((end-start) / CLOCKS_PER_SEC);

	FILE *outptasks;
	outptasks=fopen("taskres.dat", "a");
	fprintf(outptasks, "\nTask %d:\nmeanx = %.8Lf\t meanxsqu = %.8Lf\t meanspeed = %.8lf\t mu = %.8lf\t deff = %.8lf\n", taskid, 
			    tcoeff.meanx, 
			    tcoeff.meanxsqu,
			    tcoeff.meanspeed, 
			    tcoeff.mu, 
			    tcoeff.deff);

	fprintf(outptasks, "Time of simulation for task %d: %d days %dhours %dmin %dsec\n", taskid, 
                                                                                        timediff/(3600*24), 
                                                                                        (timediff/3600)%24, 
                                                                                        (timediff/60)%60, 
                                                                                        timediff%60);
	fclose(outptasks);
  }

}


void print_resallthreads(
			 long double msdall, 
                         double meanspeedall, 
			 double muall, 
                         double deffall, 
                         long double meanxall, 
                         long double meanxsquall, 
                         long double thirdcumall)
/**
 * prints results of simulation to file at the end of the simulation
 */
{
  FILE *outp;
  outp=fopen(printres.fname ,"a");
  fprintf(outp, "\n\nAverage of all Threads:\n\nmeanx = %.5Lf\t meanspeed = %.5lf\t mu = %.5lf\n\n", 
          meanxall/SimParams.numtasks, 
          meanspeedall/SimParams.numtasks, 
	  muall/SimParams.numtasks);
  fclose(outp);
           
  FILE *outpmom;
  outpmom=fopen(printres.fnamemom ,"a");
  fprintf(outpmom, "\n\nAverage of all Threads:\n\nmeanx = %.5Lf\t meanxsqu = %.5Lf\t msdall = %.5Lf\t deff = %.5lf\t thirdcum = %.5LF\n\n", 
          meanxall/SimParams.numtasks, 
	  meanxsquall/SimParams.numtasks, 
          msdall/SimParams.numtasks,
	  deffall/SimParams.numtasks,  
	  thirdcumall/SimParams.numtasks);
  fclose(outpmom);
}

void print_muoverf(double muall, double deffall, char *namefile)
{
/**
 * prints results to file outside of working directory
 */
  char fnamemu[60];

  sprintf(fnamemu, "../muoverfpos_R_%.2lf_setnumb_%d.dat", R_CONF, SimParams.setnumb);

  FILE *outmu;
  outmu=fopen(fnamemu, "a");
  fprintf(outmu, "%.3lf\t %.6lf\t %.6lf\t %s\n", SimParams.F, muall/SimParams.numtasks, deffall/SimParams.numtasks, namefile);                       
  fclose (outmu); 

}

double **calloc_2Ddouble_array(int m, int n){
/**
 * Function that allocates memory for a 2 dimensional array of doubles
 */
double **array;
int i;

  array = calloc(m, sizeof(double));
  for(i = 0; i < m; i++){
  	array[i] = calloc(n, sizeof(double));
  }

  return array;
}





int main (int argc, char **argv){
/** main function of Brownian motion simulation
 *
 */

  /*
   * initialize parameters for time measurement and measure time
   */ 
  clock_t prgstart; 
  prgstart = clock(); 
  
//  double u, v, x, y, yo, xo; 
  double yue, distx, disty, dist; 
  double f_cut;
//  double f_cut,  
  double deffall, meanspeedall, muall;
//  double muabbo, deffabbo;
//  double t, dt; 
//  double xtest, ytest;
  unsigned int xcountercheck, ycountercheck, twodcountercheck;
  int i, j;
  long double  meanxall, meanxsquall, msdall, thirdcumall;
//  long double  sqrt_flucts, f_dt;
//  int abb, abbdeff;
//  int int_ind, kset, shiftind; 
  int tasks = 1;
  int taskid = MASTER;
  int setn_per_task;
  int setn;
  char *namefile;
  char *conprfx;
  char *intprfx;
  bool ParameterFlag;
/*  bool PrintRes;
  bool TestRes;*/

  /*
   * width of bins used for spatial discretization for position histograms   
   */
  double binx;
  double biny;
  double bin2d; 
  
  /* init state for print functions to check if
   * header line in result file has been printed
   */
  printres.state = 0; 
  
  /*
   * read external force and number of interacting particles per set from command line arguments
   */ 
  printf("\nargc: %d\n\n", argc);
  if(argc > 3){
	SimParams.F = atof(argv[1])*L;
	SimParams.setnumb = atof(argv[2]);
  }
  else{
 	return -1;
  }
  
  SimParams.time_step = PARAMS_time_step(B, R_INT); 
/*  sqrt_flucts = sqrt(2*BOTTRAD*SimParams.time_step);
  f_dt = SimParams.F*SimParams.time_step;*/

  
  PARAMS_init();

  binx = 0.02*L;
  biny = 2.0*binx*MAX_HALF_WIDTH/L;
  bin2d = 0.05; 


  /*
   * get name of confinement and type of intra particle interaction
   * and provide them to function that creates working directory
   */
  conprfx = CONF_prfx();
  intprfx = prfx_int();
  namefile = CODEHAND_makedirectory(conprfx, intprfx);
  chdir(namefile);


  /**
   * calculate number of samples of interacting particles
   */ 
  setn = (int) SimParams.N/SimParams.setnumb;
  
  
  /*
   * initialize arrays where x- and y-coordinates of particles are stored in
   */ 
  printf("\ninitialize arrays for positions and interactions\n"); 
  double **positionx;
  double **positiony;
  double **fintxarray;
  double **fintyarray;
  double **xstart;

  
  positionx = calloc_2Ddouble_array(setn, SimParams.setnumb);
  positiony = calloc_2Ddouble_array(setn, SimParams.setnumb);
  
  fintxarray = calloc_2Ddouble_array(setn, SimParams.setnumb);
  fintyarray = calloc_2Ddouble_array(setn, SimParams.setnumb);
  
  xstart = calloc_2Ddouble_array(setn, SimParams.setnumb);
  
  char fnamex [60];
  char fnamey [60];
  char fname2d [60];
  char fnamespecs [60];

  CODEHAND_copy_main();
  PARAMS_copycode(); 
  CONF_copycode();
  copycode_int();
  RES_copycode();
  CODEHAND_copycode();

  sprintf(fnamespecs, "simulation_specs.dat");
  PARAMS_basic(fnamespecs);
  CONF_specs(binx, biny, bin2d);
  f_cut = intforce(INT_CUTOFF, INT_CUTOFF);
  specs_int(f_cut);


# ifdef MPI_ON
	  MPI_Init (&argc,&argv);
	  MPI_Comm_size (MPI_COMM_WORLD, &tasks);
	  MPI_Comm_rank (MPI_COMM_WORLD, &taskid); 
# endif

  SimParams.numtasks = tasks;
  printf("\n numtasks: %d\n", SimParams.numtasks);
  
  ParameterFlag = PARAMS_check_consistency();
  if(ParameterFlag == false){
	return -1;
  }	
  
  /**
   * calculate number of samples of interacting particles per task
   */ 
  setn_per_task = (int) SimParams.N/(SimParams.setnumb*SimParams.numtasks);

  if(taskid == MASTER){
  	  if(SimParams.numtasks > 1){
		  FILE *outptasks;
		  outptasks=fopen("taskres.dat", "a");
		  fprintf(outptasks, "\nFile shows results of all tasks if code is parallelized:\n\n");
		  fclose(outptasks);
          } 

  }
  
  /* Initialize pointer r as interface to gls random functions */  
  gsl_rng *r = gsl_rng_alloc (gsl_rng_mt19937);  
  /* set time(NULL + taskid as 'random' seed */ 
  gsl_rng_set(r, time(NULL) + taskid);
  
  
  /* 
   * Initialize particle positions 
   */
  SIM_init_positions(setn_per_task, positionx, positiony, xstart, r);


  /* Initialize inter-particle forces */
  SIM_init_interactions(setn_per_task, positionx, positiony, fintxarray, fintyarray);

  printf("\ntask ID:\t %d -> particle positions and forces fixed\n", taskid);
 
  /*
   * Core of Simulation: particles are propagated until criteria
   * for equilibration are met
   */  
  SIM_simulation_core(setn_per_task,
		      setn,
		      taskid, 
                      positionx, 
                      positiony, 
                      xstart, 
                      fintxarray,
                      fintyarray,
                      r);
 
  /*Merge quantities that were calculated in separated MPI threads*/ 
#  ifdef MPI_ON
	  MPI_Reduce(&tcoeff.meanx, &meanxall, 1, MPI_LONG_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
	  MPI_Reduce(&tcoeff.meanxsqu, &meanxsquall, 1, MPI_LONG_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
	  MPI_Reduce(&tcoeff.msd, &msdall, 1, MPI_LONG_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
	  MPI_Reduce(&tcoeff.thirdcum, &thirdcumall, 1, MPI_LONG_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
	  MPI_Reduce(&tcoeff.meanspeed, &meanspeedall, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
	  MPI_Reduce(&tcoeff.mu, &muall, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
	  MPI_Reduce(&tcoeff.deff, &deffall, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
#  else
        meanxall = tcoeff.meanx;
        meanxsquall = tcoeff.meanxsqu;
        msdall = tcoeff.msd;
        thirdcumall = tcoeff.thirdcum;
        meanspeedall = tcoeff.meanspeed;
        muall = tcoeff.mu;
        deffall = tcoeff.deff;
#  endif
  
  print_runtime_threads(prgstart, SimParams.numtasks, taskid);
  
  sprintf(fnamex, "meanx_Histogram_F_%.3lf.dat", SimParams.F);
  xcountercheck = RES_histogramm_mpi_reduce(setn_per_task, 
					0, 
					L, 
					binx, positionx, 
					fnamex, taskid);

  sprintf(fnamey, "meany_Histogram_F_%.3lf.dat", SimParams.F);
  ycountercheck = RES_histogramm_mpi_reduce(setn_per_task, 
                                        MAX_HALF_WIDTH, 
                                        2*MAX_HALF_WIDTH, 
                                        biny, positiony, 
                                        fnamey, 
                                        taskid);

  sprintf(fname2d, "meanpos_Histogram2d_F_%.3lf.dat", SimParams.F);
  twodcountercheck =  RES_histogramm2d_mpi_reduce(setn_per_task, 
					      bin2d, 
                                              positionx, 
					      positiony, 
                                              fname2d, 
					      taskid);


  if(taskid == MASTER){
 	   print_positions(setn_per_task, positionx, positiony);
           RES_print_countercheck(xcountercheck, ycountercheck, twodcountercheck, fnamespecs);
           print_resallthreads(msdall, meanspeedall, muall, deffall, meanxall, meanxsquall, thirdcumall);
	   print_muoverf(muall, deffall, namefile);
            
       //    CODEHAND_delerrorfiles();           
            
	   print_runtime(prgstart);
           
  }
  
  free(positionx);
  free(positiony);
  free(fintxarray);
  free(fintyarray);
  free(xstart);
  
  #ifdef MPI_ON
	  MPI_Finalize();
  #endif

  gsl_rng_free(r);
  r = NULL;

  chdir("../");
  return 0;
}
