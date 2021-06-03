/*Calculation of the averagespeed for a brownian partical with finite radius in 2d in a confinement, 
trajektories are calculatet parallel*/

#include "code_handling.h"
#include "simulation_core.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>


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
  /* effective diffusion coefficient accumulated from all threads
   * if MPI parallelization is employed.
   */ 
  double deffall;
  /* mean speed of all particles accumulated from all
   * threads if MPI is employed
   */
  double meanspeedall;
  /* non-linear mobility accumulated from all threads
   * if MPI is employed.
   */
  double muall;
  /* mean position in x-direction of all
   * particles accumulated from all threads if 
   * MPI is employed.
   */
  long double  meanxall;
  /* mean squared position <x^2> accumulated from
   * all threads if MPI is employed.
   */ 
  long double meanxsquall;
  /* mean-squared displacement <x^2> - <x>^2 accumulated
   * from all threads if MPI is employed.
   */ 
  long double msdall;
  /* third cumulant of position in x-direction 
   * accumulated from all threads if MPI is employed. 
   */
  long double thirdcumall;
  /* initialize number of MPI tasks to one to
   * get consistent behavior if MPI is not employed.
   */
  int tasks = 1;
  /* taskid of process calling main is MASTER */
  int taskid = MASTER;
  /* number of particle ensembles simulated in one task.
   * this number has to be at least one, to prevent 
   * necessary communication between the threads during 
   * simulation time.
   */
  int setn_per_task;
  /* number of 'sets' or ensembles simulated */
  int setn;
  /* pointer to variable where directory name
   * where code and results are copied to
   */
  char *namefile;
  /* prefix to mark output-files with employed confinement */
  char *conprfx;
  /* prefix to mark output-files with employed 
   * inter-particle interaction */
  char *intprfx;


  
  /* init state for print functions to check if
   * header line in result file has been printed
   */
  printres.state = 0; 
  
  /*
   * read external force and number of interacting particles per set from command line arguments
   */ 
  printf("\nargc: %d\n\n", argc);
  if(argc > 2){
	SimParams.F = atof(argv[1])*L_CONF;
	SimParams.setnumb = atof(argv[2]);
//	mpi_ind = atof(argv[3]);
  }
  else{
 	return -1;
  }
  
  SimParams.time_step = PARAMS_time_step(B, R_INT); 
	
  /*
   * initialize parameters and resulting values
   */
  PARAMS_init();
  RES_init(); 


  /*
   * get name of confinement and type of intra particle interaction
   * and provide them to function that creates working directory
   */
  conprfx = CONF_prfx();
  intprfx = INT_prfx();
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
  char fname_simparams [60];
  char fname_confparams [60];
  char fname_intparams [60];

# ifdef MPI_ON
	  MPI_Init (&argc,&argv);
	  MPI_Comm_size (MPI_COMM_WORLD, &tasks);
	  MPI_Comm_rank (MPI_COMM_WORLD, &taskid); 
# endif
  SimParams.numtasks = tasks;

  if(taskid == MASTER){
	  CODEHAND_copy_main();
	  PARAMS_copycode(); 
	  CONF_copycode();
	  INT_copycode();
	  RES_copycode();
	  CODEHAND_copycode();
	  SIM_copycode();
	  PRINT_copycode();
	  CODEHAND_copy_comp_gen_header();

	  /*filenames of simulation parameters*/
	  sprintf(fname_simparams, "parameters_simulation_overall.dat");
	  sprintf(fname_confparams, "parameters_confinement.dat");
	  sprintf(fname_intparams, "parameters_particle_interaction.dat");
	  
	  PARAMS_basic(fname_simparams);
	  CONF_specs(fname_confparams);
	  INT_specs(fname_intparams);
	  
	  printf("\n numtasks: %d\n", SimParams.numtasks);
	  
	  PARAMS_check_consistency();
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
  
  PRINT_runtime_threads(prgstart, SimParams.numtasks, taskid);
  
  sprintf(fnamex, "meanx_Histogram_F_%.3lf.dat", SimParams.F);
  histparams.xcounter = RES_histogramm_mpi_reduce(setn_per_task, 
						0, 
						L_CONF, 
						histparams.binx, positionx, 
						fnamex, taskid);

  sprintf(fnamey, "meany_Histogram_F_%.3lf.dat", SimParams.F);
  histparams.ycounter = RES_histogramm_mpi_reduce(setn_per_task, 
						MAX_HALF_WIDTH, 
						2*MAX_HALF_WIDTH, 
						histparams.biny, positiony, 
						fnamey, 
						taskid);

  sprintf(fname2d, "meanpos_Histogram2d_F_%.3lf.dat", SimParams.F);
  histparams.twodcounter =  RES_histogramm2d_mpi_reduce(setn_per_task, 
						      histparams.bin2d, 
						      positionx, 
						      positiony, 
						      fname2d, 
						      taskid);


  if(taskid == MASTER){
 	   PRINT_positions(setn_per_task, positionx, positiony);
           RES_print_countercheck(fname_confparams);
           PRINT_resallthreads(msdall, meanspeedall, muall, deffall, meanxall, meanxsquall, thirdcumall);
	   PRINT_muoverf(muall, deffall, namefile);
            
       //    CODEHAND_delerrorfiles();           
            
	   PRINT_runtime(prgstart, fname_simparams);
           
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
