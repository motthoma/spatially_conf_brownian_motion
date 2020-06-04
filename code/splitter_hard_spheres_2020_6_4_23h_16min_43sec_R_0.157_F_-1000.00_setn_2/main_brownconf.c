/*Calculation of the average speed for a brownian partical with finite radius in 2d in a confinement, 
trajektories are calculatet parallel*/


#include "par_sim.h"
#include "comp_gen_header.h"
#include "results_transport.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define MASTER 0    



/**
 * structure which contains variables used in the
 * printing functions which print the resulting
 * transport coefficients in .dat files
 */
struct PrintResults{
	/*state index to monitor if header line
         *or ongoing value has to be printed
         */
	int state;
        /*array to store name of .dat file containing mobility and other coeffs*/
	char fname [60];
        /*array to store name of .dat file containing moments of positions*/
	char fnamemom [60];

} printres;

/**
 * Function that caluclates transport coefficients such as
 * mobility, means-squared displacement or thir moment of
 * position. Function takes relative x-positions and values
 * of shift to first call function for calculation of individual
 * absolute particle positions and subsequently ensemble averages 
 * of the first three moments of the position.
 */ 
void calc_transpcoeffs(
		       int setn_per_task, 
		       double t,
		       long int **posshift,
		       long int **negshift,
		       double **posx,
		       double **x_init){

  long double totalshift = 0;
  long double abspos_part = 0;
  long double xges = 0;
  long double xgessquare = 0;
  long double xgesqub = 0;

  int j;
  int kset;
  for(j = 0; j < setn_per_task; j++){
	  for(kset = 0; kset < SimParams.setnumb; kset++){
		 totalshift = posshift[j][kset] - negshift[j][kset];

		 /*absolute position of individual particle*/
		 abspos_part = (posx[j][kset] - x_init[j][kset]) - totalshift*L;
                 /*sum of positions of all particles for later averaging*/
		 xges += abspos_part;
		 /*sum of x^2 over all particles*/
		 xgessquare += powl(abspos_part, 2);
                 /*sum of x^3 over all particles*/
		 xgesqub += powl(abspos_part, 3);
 
	  } 
  }
  
  /*calculate transport coefficients, where the averages are ensemble averages*/
  /*mean position <x> of all particles*/
  tcoeff.meanx = xges/SimParams.N*SimParams.numtasks; 
  /*mean square position <x^2> of all particles*/
  tcoeff.meanxsqu = xgessquare/SimParams.N*SimParams.numtasks; 
  /*mean squared displacement <x^2> - <x>^2 */
  tcoeff.msd = tcoeff.meanxsqu - tcoeff.meanx*tcoeff.meanx;
  /*means speed <v> = <x>/t */
  tcoeff.meanspeed = tcoeff.meanx/t;  
  /*mobility mu = <v>/F */
  tcoeff.mu = tcoeff.meanspeed/SimParams.F;
  /*effective diffusion D_eff = (<x^2>-<x>^2)/(2*t*b/r) */
  tcoeff.deff = tcoeff.msd/(2*t*BOTTRAD); 
  /*third moment of position <x^3> */
  tcoeff.meanxqub = xgesqub/SimParams.N*SimParams.numtasks;
  /*thrid cumulant of poisition <x^3> - 3<x><x^2> + 2<x>^3 */
  tcoeff.thirdcum = tcoeff.meanxqub - 3*tcoeff.meanx*tcoeff.meanxsqu + 2*powl(tcoeff.meanx, 3);
}


char* makedirectory(char *confprfx, char *intprfx){             
/**
 * creates directory, where code and data is transferred to
 * moves to this directory as working directory
 */
  time_t tnow;
  struct tm *tmnow;
  time(&tnow);
  tmnow = localtime(&tnow);

  char *e = malloc(100*sizeof(char));

  sprintf(e, "%s_%s_%d_%d_%d_%dh_%dmin_%dsec_R_%.3lf_F_%.2lf_setn_%d", confprfx, intprfx, tmnow->tm_year + 1900, tmnow->tm_mon+1, tmnow->tm_mday, tmnow->tm_hour,tmnow->tm_min,tmnow->tm_sec, R_CONF, SimParams.F, SimParams.setnumb);
  /**
   * 0700 is modus for determining access rights to dir
   */
  mkdir(e, 0700);
  chdir(e);
  return e;
}

void copy_main(){
/**
 * copies main_brownianconf.c to directory created by 'makedirectory'
 */

char copycode[200];

  sprintf(copycode, "cp ../main_brownconf.c ../masterinteract.py ./");
  system(copycode);

}

void delerrorfiles(){
/**
 * deletes error and log files used on albeniz
 */
char delerrorfile[100];
  sprintf(delerrorfile, "rm ../F_%.2lf_sn_%.0d_clintparallel* ",SimParams.F, SimParams.setnumb);
  system(delerrorfile);

}

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

void print_results_over_time(double t, 
                             int abb, 
                             int abbdeff)
{
       /**
	* function to print online results of simulation over time 
	*/
        if(printres.state == 0)	{    
		sprintf(printres.fname, "muovert_F_%.3lf.dat", SimParams.F);
		FILE *outp;
		outp = fopen(printres.fname ,"w");
		fprintf(outp, "#time\t meanx: <x>\t meanspeed: <v>\t mu\t abb\t abbdeff\n");
		fclose(outp);

		sprintf(printres.fnamemom, "momsovert_F_%.3lf.dat", SimParams.F);
		FILE *outpmom;
		outpmom=fopen(printres.fnamemom ,"w");
		fprintf(outpmom, "#time\t meanx: <x>\t Meansqdist: <x^2> - <x>^2\t  deff: (<x^2> - <x>^2)/(2t)\t third cumulant: <x^3>-3<x^2><x>+2<x>^3\n");
		fclose(outpmom);

		printres.state = 1;
	}
	else{
		FILE *outp;
		outp=fopen(printres.fname ,"a");
		fprintf(outp, "%.6f\t %.4Lf\t %.5lf\t %.4lf\t  %d\t %d\n", t, tcoeff.meanx, tcoeff.meanspeed, tcoeff.mu, abb, abbdeff);
		fclose(outp);

		FILE *outpmom;
		outpmom=fopen(printres.fnamemom ,"a");
		fprintf(outpmom, "%.6f\t %.6Lf\t %.5Lf\t %.5Lf\t %.5lf\t %.5Lf\n", t, tcoeff.meanx, tcoeff.meanxsqu, tcoeff.msd, tcoeff.deff, tcoeff.thirdcum);
		fclose(outpmom);
		
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

int histogramm_mpi_reduce(int m, 
                          double backshift, 
                          double length, 
                          double bin, 
                          double **positions, 
                          char *fname, 
                          int taskid)
{
/**
 * Stores 1 dimensional histogram in file named fname.
 * During a spatial sweep through all length/bin slices from position 
 * 'backshift' to position 'length', all particle positions are checked 
 * and the number of particles in the respective slice is increased if 
 * a particle is detected. The number counter of all particles in the 
 * slice is stored in the file together with the upper and lower 
 * boundary of the slice.
 */

  int i;
  int j;
  int k;
  int bin_n; 
  int countercheck = 0;
  int counter = 0;
  int counterall = 0;
  
  FILE *outp;
  outp = fopen(fname, "w");
  fclose(outp);      

  bin_n = (int) (length/bin);
  for(i = 0; i < bin_n; i++){ 
	  counter = 0;
	  counterall = 0;
        
	  for(j = 0; j < m; j++){
		for(k = 0; k < SimParams.setnumb; k++){
		   if(((i*bin - backshift) <= positions[j][k]) && (positions[j][k] <= (i+1)*bin - backshift)){
			counter++;
		   }

		}
	  }
	  #ifdef MPI_ON
		  MPI_Reduce(&counter, &counterall, 1, MPI_INTEGER, MPI_SUM, MASTER, MPI_COMM_WORLD);
	  #else
		  counterall = counter;
	  #endif
	  if(taskid == MASTER){
                  countercheck += counterall;
		   
		  outp = fopen(fname, "a");
		  fprintf(outp,"%f\t %f\t %d\t %f\n", i*bin - backshift, (i+1)*bin - backshift, counterall, counterall/(SimParams.N*bin));
		  fclose(outp);      
	  }
  }

  return countercheck;

}

int histogramm2d_mpi_reduce(int m, 
                            double bin2d, 
                            double **positionsx, 
                            double **positionsy,  
                            char *fname, 
                            int taskid)
{
/**
 * Stores 2 dimensional histogram in file named fname.
 * During a spatial scan over all nx*ny rectangular fields, all particle
 * positions are checked and the number of particles in the respective
 * field is increased if a particle is detected. The number counter
 * of all particles in the field is stored in the file together
 * with the upper and lower boundaries in both directions of the field.
 */

  int i;
  int j;
  int hx;
  int hy;
  
  int bin_nx = (int) (L/bin2d);
  int bin_ny = (int) (2*MAX_HALF_WIDTH/bin2d);
  int twodcountercheck = 0; 
  int twodcounter = 0;   
  int twodcounterall = 0;   
  
  FILE *outp;
  outp = fopen(fname, "w");
  fclose(outp);      
 
   for(hx = 0; hx <= bin_nx; hx++){ 
	  for(hy = 0; hy <= bin_ny; hy++){   
                  twodcounter = 0;   
                  twodcounterall = 0;   
		  for(i = 0; i < m; i++){
			for(j = 0; j < SimParams.setnumb; j++){
			   if((hx*bin2d <= positionsx[i][j]) && (positionsx[i][j] <= (hx+1)*bin2d)){
				   if(((hy*bin2d - MAX_HALF_WIDTH) <= positionsy[i][j]) && (positionsy[i][j] <= (hy+1)*bin2d - MAX_HALF_WIDTH)){
					   twodcounter++;
				   }
			   }   

			}
		  }

		  #ifdef MPI_ON
			  MPI_Reduce(&twodcounter, &twodcounterall, 1, MPI_INTEGER, MPI_SUM, MASTER, MPI_COMM_WORLD);
		  #else		  
			  twodcounterall = twodcounter;
		  #endif
		  if(taskid == MASTER){
			  twodcountercheck += twodcounterall;

			  outp = fopen(fname, "a");
			  fprintf (outp, "%f\t%f\t%f\t%f\t\t%d\t\t%f\n", hx*bin2d, (hx+1)*bin2d, hy*bin2d - MAX_HALF_WIDTH, (hy+1)*bin2d - MAX_HALF_WIDTH, twodcounterall, twodcounterall/(SimParams.N*bin2d*bin2d));
			  fclose(outp);
		  }
	  }
  }
  return twodcountercheck;
}
  
void init_particle_pos(int setn_per_task, 
                       double **positionx, 
                       double **positiony, 
                       double **xstart, 
                       gsl_rng *r)
{
/**
 * function that initializes the particle positions. The particles are uniformly distributed
 * over the whole period length L in a strip of width SimParams.initwidth or the width of
 * the channel if it is smaller. For interacting particles with hard core radius R_int, 
 * configurations which involve an overlap of two particles are rejected. 
 */
  int j;
  int kset;
  int ktest;
  double xo;
  double yo;
  double yue;
  double distx;
  double disty;
  double distinit;
  bool PosValidInit;

  PosValidInit = false;
  for(j = 0; j < setn_per_task; j++){
	  for(kset = 0; kset < SimParams.setnumb; kset++){
		  do{
			  positionx[j][kset] = gsl_rng_uniform(r)*L;
			  positiony[j][kset] = (2*gsl_rng_uniform(r) - 1)*SimParams.initwidth;
			    
			  xo = positionx[j][kset];
			  yo = positiony[j][kset];
		 
			  yue = yuef_ext(xo, yo);
			 
			  PosValidInit = true; 
			  if(fabs(positiony[j][kset]) >= yue) PosValidInit = false;
			  if((kset > 0) && (PosValidInit == true)){
					    
				  for(ktest = 0; ktest < kset; ktest++){
					  distx = xo - positionx[j][ktest];			
					  disty = yo - positiony[j][ktest];			
					  distinit = sqrt(distx*distx + disty*disty);
					  if(distinit <= 2*R_INT) PosValidInit = false;
				 }
			 }  
			  
		  }while(PosValidInit == false);
		
		  xstart[j][kset] = xo;
	  }
  }	
}

void init_particle_int(int setn_per_task, 
                       double **positionx, 
                       double **positiony,
                       double **fintxarray,
                       double **fintyarray)
{
/**
 * Function that calculates depending on the positions of the particles
 * the initial particle-particle interaction force if such a force e.g.
 * given by a Lennard-Jones interaction is given.
 */
  int j;
  int kset;
  int ktest;
  double distx;
  double disty;
  double dist;
  double fintxpair;
  double fintypair;
  double fintx;
  double finty;
                       
  for (j = 0; j < setn_per_task; j++){  
          for(kset = 0; kset < SimParams.setnumb; kset++){
	          fintx = 0;
		  finty = 0;
                  for(ktest = 0; ktest < SimParams.setnumb; ktest++){
                          if(kset != ktest){
                                  distx = positionx[j][kset] - positionx[j][ktest];
                                  disty = positiony[j][kset] - positiony[j][ktest];
                                  dist = sqrt(distx*distx + disty*disty);
                                  if(dist <= INT_CUTOFF){
                                          fintxpair = intforce(distx, dist);
                                          fintypair = intforce(disty, dist);
                                          fintx += fintxpair;
                                          finty += fintypair;
                                  }
                          }
                  }
                  fintxarray[j][kset] = fintx;
                  fintyarray[j][kset] = finty;
          }
  }
}

void print_hist_countercheck(int xcheck, int ycheck, int twodcheck, char *fname_specs)
{
/**
 * Function that prints the result of the counter checks available from the 
 * histogram functions. There, the total number of particles is counted 
 * in order to check if all particles are represented in the histogram
 * and situated in the confinement.
 */
	FILE *outp;
	outp=fopen(fname_specs, "a");
	fprintf(outp, "\n\nxcountercheck: %d\nycountercheck: %d\ntwodcountercheck: %d\n\n", xcheck, ycheck, twodcheck);
	fclose(outp);

	if((xcheck != SimParams.N) || (ycheck != SimParams.N) || (twodcheck != SimParams.N)){
	   printf("Error in Histogrammcounter!\n");
	}
}

double reset_pos_time(int setn_per_task, long int **posshift, long int **negshift) 
{
/**
 * Function to reset the simulated system time t and the positions
 * to originate values in order to skip transient effects from
 * start of simulation.
 */
	int i, j;
	double t = 0;
	for(i = 0; i < setn_per_task; i++){
	  for(j = 0; j < SimParams.setnumb; j++){
		 posshift[i][j] = 0; 
		 negshift[i][j] = 0;

	  } 
	}
	return t;
}


void adapt_posshifts(int shiftind, int i, int j, long int **posshift, long int **negshift)
{
/**
 * Function to update shifts which are monitored to calculate 
 * absolute position in x-direction from simulation with cyclic
 * boundary conditions
 */
	/*
	 *shift particle in positive direction if 
	 *position is 'left' of considered channel period   
	 */
	if(shiftind < 0){
	  posshift[i][j]++; 
	} 
	/*
	 *shift particle in negative direction if 
	 *position is 'right' of considered channel period   
	*/
	if(shiftind > 0){
	  negshift[i][j]++;
	}
}

int update_equcounter(double tran_quant, double tran_quanto, double accurarcy, int equcounter)
{
/**
 * Function for the update of the counter that is used to monitore the
 * equilibration of the system. A transport quantity tran_quant is
 * compared with a previous value. If the difference is  below the    
 * demanded accurarcy, the value of the counter is increased.
 * If the difference is larger than twice the accurarcy the counter
 * is decreased
*/ 
	if(fabs(tran_quant - tran_quanto) <= accurarcy){
	  equcounter++;
	}

	if(fabs(tran_quant - tran_quanto) >= accurarcy){
	  equcounter--;
	  if(equcounter < 0) equcounter = 0;
	}
	return equcounter;
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

long int **calloc_2Dlint_array(int m, int n){
/**
 * Function that allocates memory for a 2 dimensional array of long ints
 */
long int **array;
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
  
  double u, v, x, y, yo, xo, yue, distx, disty, dist, f_cut,  deffall, meanspeedall, muall;
  double muabbo, deffabbo;
  double t, dt, xtest, ytest;
  unsigned int xcountercheck, ycountercheck, twodcountercheck;
  int i, j;
  long double  meanxall, meanxsquall, msdall, thirdcumall;
  long double  sqrt_flucts, f_dt;
  int abb, abbdeff, int_ind, kset, shiftind; 
  int tasks = 1;
  int taskid = MASTER;
  int setn_per_task;
  int setn;
  char *namefile;
  char *conprfx;
  char *intprfx;
  bool ParameterFlag;
  bool PosValid;
  bool PrintRes;
  bool TestRes;

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
	SimParams.F = atof(argv[1]);
	SimParams.setnumb = atof(argv[2]);
  }
  else{
 	return -1;
  }
  
  SimParams.time_step = PARAMS_time_step(B, R_INT); 
  sqrt_flucts = sqrt(2*BOTTRAD*SimParams.time_step);
  f_dt = SimParams.F*SimParams.time_step;

  
  PARAMS_init();
  TCOEF_init();

  binx = 0.02*L;
  biny = 2.0*binx*MAX_HALF_WIDTH/L;
  bin2d = 0.05; 


  /*
   * get name of confinement and type of intra particle interaction
   * and provide them to function that creates working directory
   */
  conprfx = prfx_conf();
  intprfx = prfx_int();
  namefile = makedirectory(conprfx, intprfx);
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

  printf("arrays for positions and interactions are initialized!\n"); 


  long int **negshift;
  long int **posshift;
  
  negshift = calloc_2Dlint_array(setn, SimParams.setnumb);
  posshift = calloc_2Dlint_array(setn, SimParams.setnumb);

  printf("arrays for negshift and posshift are initialized!\n"); 
  
  char fnamex [60];
  char fnamey [60];
  char fname2d [60];
  char fnamespecs [60];

  copy_main();
  PARAMS_copycode(); 
  copycode_conf();
  copycode_int();
  
  sprintf(fnamespecs, "simulation_specs.dat");
  PARAMS_basic(fnamespecs);
  specs_conf(binx, biny, bin2d);
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
  
  printf("start to init positions!\n"); 
  
  /* 
   * Initialize particle positions 
   */
  init_particle_pos(setn_per_task, positionx, positiony, xstart, r);

  printf("positions initialized!\n"); 

  /* Initialize inter-particle forces */
  init_particle_int(setn_per_task, positionx, positiony, fintxarray, fintyarray);

  printf("\ntask ID:\t %d -> particle positions and forces fixed\n", taskid);
  //printf("\nepsilon:\t %lf\n", EPS_L);
   
 
  /* Perform simulation steps until equilibration is reached */
  t = 0;
  dt = SimParams.time_step;  
  i = 1;      
  abb = 0;	
  abbdeff = 0;	
  muabbo = 0;	
  deffabbo = 0;	
  /*loop over simulation steps.
   * loop is stopped when criterion for equilibration is fulfilled.
   */
  do{	
	
	  t += dt;
	  i++;
	  /* Loop over trajectories */
	  for (j = 0; j < setn_per_task; j++){	
		  for(kset = 0; kset < SimParams.setnumb; kset++){


			  xo = positionx[j][kset];
			  yo = positiony[j][kset];
			  
			  fintx = fintxarray[j][kset];
			  finty = fintyarray[j][kset];
		  
		 	  /* 
 			   * Perform simulation steps until valid step where no particles overlap
 			   * and particles are within channel is obtained 
 			   */
			  do{
                                 /*
                                  *Create random numbers for x- and y-component of noise
                                  */  
				  u = gsl_ran_gaussian_ziggurat(r,1);
				  v = gsl_ran_gaussian_ziggurat(r,1);  
				 /* 
                                  *Update x- and y-component of position according to
                                  *stochastic Euler for Langevin equation
                                  */                       
				  x = xo + f_dt + fintx*dt + sqrt_flucts*u;
				  y = yo + finty*dt + sqrt_flucts*v;
			
                                  /*
                                   *Calculate value of confinement boundary at current
                                   *position x (y-value is needed for non-analytic treatment
                                   *of channels with cosine shape
                                   */  
				  yue = yuef_ext(x,y);
				    
				  PosValid = true;
				  /*Check if particle is within effective boundary*/  
				  if (fabs(y) > yue) PosValid = false;
				  /*Check if bottleneck is passed correctly*/	
				  if ((x < 0) && ((fabs(yo) >= B-R_CONF) || (fabs(y) >= B-R_CONF))) PosValid = false;
				  if ((x > L) && ((fabs(yo) >= B-R_CONF) || (fabs(y) >= B-R_CONF))) PosValid = false;

				  if (PosValid == true){			
					  shiftind = 0;
					  if(x < 0){
						   shiftind = -1;
						   x += L;
					  }
					  if(x > L){
						   shiftind = 1;
						   x -= L;
					  }
					  /*simulate particle-particle interaction*/ 
 					  if (SimParams.setnumb > 1){ 
						    fintx = 0;
						    finty = 0;
 						    for(int_ind = 0; int_ind < SimParams.setnumb; int_ind++){	
							 
							  xtest = positionx[j][int_ind]; 
							  ytest = positiony[j][int_ind];
	 
							  if(int_ind != kset){
								  distx = x - xtest;
								  disty = y - ytest;

								  /* 
                                                                   * search relevant distance according
 								   * to minimum image conversion 
 								   */
								  if (abs(distx) > 0.5*L){ 
									distx = distx - L*(distx/abs(distx));
								  }
								  dist = sqrt(distx*distx + disty*disty);
								  
							          if(dist <= 2*R_INT){ 
								  	break;
								  }
								  
								  if(dist <= INT_CUTOFF){
									  fintxpair = intforce(distx, dist);
									  fintypair = intforce(disty, dist);
									  fintx += fintxpair;
									  finty += fintypair;

								  }
							  }

						    /*close loop over particles of interacting ensemble*/
						    } 
					    }		    
				  }
                              
			  }while(PosValid == false);
			  
		          if(shiftind != 0){
				  adapt_posshifts(shiftind, j, kset, posshift, negshift);
                          }
		  
                         /*
                          * Update arrays with positions and inter particle forces
                          */     
		          positionx[j][kset] = x; 
		          positiony[j][kset] = y;
			  fintxarray[j][kset] = fintx;
			  fintyarray[j][kset] = finty;

		}
 	  /*
           * Close loop over trajectories 
           */
	  }
	
       	  /*
	   *Initialize reference values of mobility and diffusion coefficients
	   *for later judgement of equilibration process.
	  */
	  if((i > SimParams.stepnumb - SimParams.testab) && (i <= SimParams.stepnumb - SimParams.testab + 1)){ 
		  muabbo = tcoeff.mu;
		  deffabbo = tcoeff.deff;
	  }
	
          /*
           * Test progress of equilibration and plot results at certain 
           * simulation steps i
           */ 
          PrintRes = false;
	  if(i % SimParams.plotpoints == 0){
	  	PrintRes = true;
	  }
	  TestRes = false;
	  if((i > SimParams.stepnumb) && (i % SimParams.testab == 0)){
     		TestRes = true;
          }

	  if((TestRes == true) || PrintRes == true){ 

		  /*
		   * call function for calculation of transport coefficients 
		   * such as mobility or mean-squared displacement 
		   * */ 
		  calc_transpcoeffs(setn_per_task, 
				    t, 
				    posshift,
				    negshift,
				    positionx,
				    xstart);
	  
			 
                  /*
                   * Update of the equilibration counter that are used to monitore the
                   * equilibration of the mobility and diffusivity
                   */  
		  if(TestRes == true){ 
 			  abb = update_equcounter(tcoeff.mu, muabbo, SimParams.accur, abb);

			  if(tcoeff.deff > 1.0){  
				  
 			  	  abbdeff = update_equcounter(tcoeff.deff, deffabbo, SimParams.deffaccur*deffabbo, abbdeff);
			  }

			  else{ 
 			  	  abbdeff = update_equcounter(tcoeff.deff, deffabbo, SimParams.deffaccur, abbdeff);
			  }		
			  muabbo = tcoeff.mu;    
			  deffabbo = tcoeff.deff;
                  }
 		  /*
                   * Plot results to check progress of equilibration 
                   */
		  if((PrintRes == true) && (taskid == MASTER)){
			  print_results_over_time(t, 
						  abb, 
						  abbdeff);

		  }
	  }

	  /*
           *  Reset position and time information to truncate transient effects from small times 
           */
	  if(i == SimParams.reset_stepnumb){
		  t = reset_pos_time(setn_per_task, posshift, negshift); 
	  } 
 
  /* 
   * Closes while loop over simulation steps if criteria for equilibration are fulfilled
   */
  }while((abb < SimParams.numbtest) || (abbdeff < SimParams.numbtest));
 
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
  xcountercheck = histogramm_mpi_reduce(setn_per_task, 
					0, 
					L, 
                                        binx, positionx, 
                                        fnamex, taskid);

  sprintf(fnamey, "meany_Histogram_F_%.3lf.dat", SimParams.F);
  ycountercheck = histogramm_mpi_reduce(setn_per_task, 
                                        MAX_HALF_WIDTH, 
                                        2*MAX_HALF_WIDTH, 
                                        biny, positiony, 
                                        fnamey, 
                                        taskid);

  sprintf(fname2d, "meanpos_Histogram2d_F_%.3lf.dat", SimParams.F);
  twodcountercheck =  histogramm2d_mpi_reduce(setn_per_task, 
					      bin2d, 
                                              positionx, 
					      positiony, 
                                              fname2d, 
					      taskid);


  if(taskid == MASTER){
 	   print_positions(setn_per_task, positionx, positiony);
           print_hist_countercheck(xcountercheck, ycountercheck, twodcountercheck, fnamespecs);
           print_resallthreads(msdall, meanspeedall, muall, deffall, meanxall, meanxsquall, thirdcumall);
	   print_muoverf(muall, deffall, namefile);
            
       //    delerrorfiles();           
            
	   print_runtime(prgstart);
           
  }
  
  free(positionx);
  free(positiony);
  free(fintxarray);
  free(fintyarray);
  free(xstart);
  free(negshift);
  free(posshift);
  #ifdef MPI_ON
	  MPI_Finalize();
  #endif

  gsl_rng_free(r);
  r = NULL;

  chdir("../");
  return 0;
}
