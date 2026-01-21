/*Calculation of the averagespeed for a brownian partical with finite radius in 2d in a confinement, 
trajektories are calculatet parallel*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/stat.h>
#include <unistd.h>

#include "array_utils.h"
#include "code_handling.h"
#include "sim_config.h"
#include "simulation_core.h"
#include "random_numb_gen.h"
#include "print_routines.h"
#include "comp_gen_header.h"

void MAIN_copy_main_code() {
    CODEHAND_copy_file_to_dest("main_brownconf.c");
    CODEHAND_copy_file_to_dest("makefile");
    CODEHAND_copy_file_to_dest("comp_gen_header.h");
}

void main_copy_ext_files(){
	  MAIN_copy_main_code();
      SIMCONFIG_copy_code();
      CODEHAND_copycode();
	  RES_copycode();
	  SIM_copycode();
	  PRINT_copycode();
	  CONF_copycode();
	  INT_copycode();
      RNG_copycode();
      UTILS_copy_code();
}

int main (int argc, char **argv){
/** main function of Brownian motion simulation */

  char *config_file_path;
  if (argc > 1) {
      config_file_path = argv[1];
  } else {
      config_file_path = "sim_params.conf";
  }

  /* initialize parameters for time measurement and measure time
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
  long double meanxall;

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

  /* flag to check if number of particles, ensemble
   * number of interacting particles and number 
   * of threats are consistent.
   */
  bool ConsistencyFlag;

  /* pointer to variable where directory name
   * where code and results are copied to
   */
  char *namefile;

  /* prefix to mark output-files with employed confinement */
  char *conprfx;
  /* prefix to mark output-files with employed 
   * inter-particle interaction
   */
  char *intprfx;

  T_EnsembleState EnsembleState;

# ifdef MPI_ON
	  MPI_Init (&argc,&argv);
	  MPI_Comm_size (MPI_COMM_WORLD, &tasks);
	  MPI_Comm_rank (MPI_COMM_WORLD, &taskid); 
# endif
  /*
   * initialize parameters and resulting values
   */
  SIMCONFIG_init(tasks, config_file_path);
 
  /*init transport parameters like mobility this simulation adresses*/ 
  RES_init(); 

  /*
   * get name of confinement and type of intra particle interaction
   * and provide them to function that creates working directory
   */
  conprfx = CONF_prfx();
  intprfx = INT_prfx();
  CODEHAND_makedirectory(conprfx, intprfx);

  /*
   * initialize arrays where x- and y-coordinates of particles are stored in
   * as well as interaction forces
   */ 
  printf("\ninitialize arrays for positions and interactions stored in enseble state\n"); 
  EnsembleState = SIM_alloc_ensemble_state(&SimParams);
  
  if(taskid == MASTER){
      main_copy_ext_files();

	  /*filenames of simulation parameters*/
	  SIMCONFIG_write_specs();
	  CONF_specs();
	  INT_specs();
	  printf("\n numtasks: %d\n", SimParams.numtasks);
	  ConsistencyFlag = SIMCONFIG_check_consistency();
	  if(ConsistencyFlag == false){
	  	return -1;
	  }

	  if(SimParams.numtasks > 1){
		  FILE *outptasks;
		  outptasks=fopen("taskres.dat", "a");
		  fprintf(outptasks, "\nFile shows results of all tasks if code is parallelized:\n\n");
		  fclose(outptasks);
          } 
  }

  /* Initialize pointer interface to gls random functions */  
  RNG_init_rng(taskid, time(NULL) + taskid);
  /* 
   * Initialize particle positions 
   */
  SIM_init_positions(&SimParams, &EnsembleState);
  /*SIM_read_in_positions(&SimParams, &EnsembleState, SimParams.setn_per_task, 
                          positionx, 
                          positiony, 
                          xstart);*/

  /* Initialize inter-particle forces */
  SIM_init_interactions(&SimParams, &EnsembleState);
  printf("\ntask ID:\t %d -> particle positions and forces fixed\n", taskid);
 
  /*
   * Core of Simulation: particles are propagated until criteria
   * for equilibration are met
   */  
  SIM_simulation_core(&SimParams, &EnsembleState, taskid); 
 
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
    meanxsquall = tcoeff.meanxsqu; msdall = tcoeff.msd;
    thirdcumall = tcoeff.thirdcum;
    meanspeedall = tcoeff.meanspeed;
    muall = tcoeff.mu;
    deffall = tcoeff.deff;
#  endif
  
  PRINT_runtime_threads(prgstart, SimParams.numtasks, taskid);
  
  char fnamex [60];
  snprintf(fnamex,
           sizeof fnamex,
           "%s/meanx_Histogram_F_%.3lf.dat",
           DestPaths.fullpath,
           SimParams.F);
  histparams.xcounter = RES_histogramm_mpi_reduce(0, 
                                                  L_CONF, 
                                                  histparams.binx,
                                                  EnsembleState.positionx, 
                                                  fnamex,
                                                  taskid);

  char fnamey [60];
  snprintf(fnamey,
           sizeof fnamey,
           "%s/meany_Histogram_F_%.3lf.dat",
           DestPaths.fullpath,
           SimParams.F);
  histparams.ycounter = RES_histogramm_mpi_reduce(MAX_HALF_WIDTH, 
                                                  2*MAX_HALF_WIDTH, 
                                                  histparams.biny,
                                                  EnsembleState.positiony, 
                                                  fnamey, 
                                                  taskid);

  char fname2d [60];
  snprintf(fname2d,
           sizeof fname2d,
           "%s/meanpos_Histogram2d_F_%.3lf.dat",
           DestPaths.fullpath,
           SimParams.F);
  histparams.twodcounter =  RES_histogramm2d_mpi_reduce(histparams.bin2d, 
                                                        EnsembleState.positionx, 
                                                        EnsembleState.positiony, 
                                                        fname2d, 
                                                        taskid);


  if(taskid == MASTER){
 	   PRINT_positions(EnsembleState.positionx,
                       EnsembleState.positiony);

       RES_print_countercheck();

       PRINT_resallthreads(msdall,
                           meanspeedall,
                           muall,
                           deffall,
                           meanxall,
                           meanxsquall,
                           thirdcumall);

	   PRINT_muoverf(muall, deffall);

	   PRINT_runtime(prgstart);
           
  }
  
  free(EnsembleState.positionx);
  free(EnsembleState.positiony);
  free(EnsembleState.fintxarray);
  free(EnsembleState.fintyarray);
  
  #ifdef MPI_ON
	  MPI_Finalize();
  #endif

  RNG_free_rng();
  return 0;
}

