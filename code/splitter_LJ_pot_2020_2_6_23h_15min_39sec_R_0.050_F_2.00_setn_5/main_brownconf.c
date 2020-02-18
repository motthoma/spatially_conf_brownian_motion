/*Calculation of the average speed for a brownian partical with finite radius in 2d in a confinement, 
trajektories are calculatet parallel*/


#include "par_sim.h"
#include "comp_gen_header.h"

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
 * structure which contains transport coefficients
 * that are calculated within individual threads.
 **/
struct TransportCoeffs{
	
        /*average of particle x-positions <x>*/
	long double meanx; 
	/*average of particle speed in x-directions <v>*/
	double meanspd;
	/*<x^2>*/
	long double meanxsqu; 
	/*<x^3>*/
        long double meanxqub;
        /*mean-squared displacement of x-position (2nd moment): <x^2> - <x>^2*/
        long double msd;
	/*third moment of x-position: <x^3> - 3*<x>*<x^2> + 2*<x>^3*/
	long double mthree;
        /*effective diffusion coeff.: (<x^2> - <x>^2)/(2*t*B/R)*/
	double deff; 
	/*non-linear mobility mu = <v>/F*/
	double mu;

} tcoeff;

/**
 * Function that caluclates transport coefficients such as
 * mobility, means-squared displacement or thir moment of
 * position. Function takes relative x-positions and values
 * of shift to first call function for calculation of individual
 * absolute particle positions and subsequently ensemble averages 
 * of the first three moments of the position.
 */ 
void calc_transpcoeffs(int numtasks,
		       int setn_per_task, 
		       double t,
		       double F, 
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
	  for(kset = 0; kset < simparams.setnumb; kset++){
		 totalshift = posshift[j][kset] - negshift[j][kset];

		 /*calculate absolute position of individual particle*/
		 abspos_part = (posx[j][kset] - x_init[j][kset])/L - totalshift;
		 xges += abspos_part;
		 xgessquare += powl(abspos_part, 2);
		 xgesqub += powl(abspos_part, 3);
 
	  } 
  }
  

  tcoeff.meanx = xges/simparams.N*numtasks; 
  tcoeff.meanxsqu = xgessquare/simparams.N*numtasks; 
  tcoeff.msd = tcoeff.meanxsqu - tcoeff.meanx*tcoeff.meanx;
  tcoeff.meanspd = tcoeff.meanx/t;  
  tcoeff.mu = tcoeff.meanspd/F;
  tcoeff.deff = tcoeff.msd/(2*t*BOTTRAD); 
  tcoeff.meanxqub = xgesqub/simparams.N*numtasks;
  tcoeff.mthree = tcoeff.meanxqub - 3*tcoeff.meanx*tcoeff.meanxsqu + 2*powl(tcoeff.meanx, 3);
}


char* makedirectory(double a, int b, char* c, char *d){             
/**
 * creates directory, where code and data is transferred to
 * moves to this directory as working directory
 */
  time_t tnow;
  struct tm *tmnow;
  time(&tnow);
  tmnow = localtime(&tnow);

  char *e = malloc(100*sizeof(char));

  sprintf(e, "%s_%s_%d_%d_%d_%dh_%dmin_%dsec_R_%.3lf_F_%.2lf_setn_%d", c, d, tmnow->tm_year + 1900, tmnow->tm_mon+1, tmnow->tm_mday, tmnow->tm_hour,tmnow->tm_min,tmnow->tm_sec, R_CONF, a, b);
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

void delerrorfiles(double a, int b){
/**
 * deletes error and log files
 */
char delerrorfile[100];

  sprintf(delerrorfile, "rm ../F_%.2lf_sn_%.0d_clintparallel* ",a, b);
  system(delerrorfile);

}

void print_runtime(clock_t start, int numtasks)
{
/**
 * prints run time of code into file
 */
 
  int timediff;
  int timediff_all;
  clock_t end = clock();
  
  timediff = (int)((end-start) / CLOCKS_PER_SEC);
  timediff_all = (int)(numtasks*(end-start) / CLOCKS_PER_SEC);

  FILE *outpspecs;
  outpspecs=fopen("muovert_specs.dat", "a");
  fprintf(outpspecs, "\nSimulationtime: %d days %dhours %dmin %dsec\n", timediff/(3600*24), (timediff/3600)%24, 
                                                                        (timediff/60)%60, timediff%60);

  fprintf(outpspecs, "Total computing time of all threads: %d days %dhours %dmin %dsec\n\n", timediff_all/(3600*24), 
                                                                                             (timediff_all/3600)%24, 
                                                                                             (timediff_all/60)%60, 
                                                                                              timediff_all%60);
  fclose(outpspecs);

}

void print_positions(int m, int n, double **posx, double **posy){
/**
 * prints particle positions to file
 */
  int i;
  int j;
  FILE *outpos;
  outpos=fopen("positions.dat", "a");
  fprintf(outpos, "#xpositions\t ypositions\n");

  for(i = 0; i < m; i++){
	for(j = 0; j < n; j++){
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
	fprintf(outptasks, "\nTask %d:\ntcoeff.meanx = %.8Lf\t tcoeff.meanxsqu = %.8Lf\t tcoeff.meanspd = %.8lf\t tcoeff.mu = %.8lf\t tcoeff.deff = %.8lf\n", taskid, 
			    tcoeff.meanx, 
			    tcoeff.meanxsqu,
			    tcoeff.meanspd, 
			    tcoeff.mu, 
			    tcoeff.deff);

	fprintf(outptasks, "Sitcoeff.mulationtime for task %d: %d days %dhours %dmin %dsec\n", taskid, 
                                                                                        timediff/(3600*24), 
                                                                                        (timediff/3600)%24, 
                                                                                        (timediff/60)%60, 
                                                                                        timediff%60);
	fclose(outptasks);
  }

}


void print_resallthreads(long double msdall, 
                         double meanspdall, 
			 double muall, 
                         double deffall, 
                         long double meanxall, 
                         long double meanxsquall, 
                         long double mthreeall,
                         int numtasks,
                         char fname[],
                         char fnamemom[])
/**
 * prints results of sitcoeff.mulation to file
 */
{
  FILE *outp;
  outp=fopen(fname ,"a");
  fprintf(outp, "\n\nAverage of all Threads:\n\nmsd = %.5Lf\t meanspd = %.5lf\t mu = %.5lf\t deff = %.5lf\n\n", 
          tcoeff.msd/numtasks,
          tcoeff.meanspd/numtasks, 
	  tcoeff.mu/numtasks, 
	  tcoeff.deff/numtasks);
  fclose(outp);
           
  FILE *outpmom;
  outpmom=fopen(fnamemom ,"a");
  fprintf(outpmom, "\n\nAverage of all Threads:\n\nmeanx = %.5Lf\t meanxsqu = %.5Lf\t mthree = %.5LF\n\n", 
          tcoeff.meanx/numtasks, 
	  tcoeff.meanxsqu/numtasks, 
	  tcoeff.mthree/numtasks);
  fclose(outpmom);
}

void print_muoverf(double F,  int numtasks, double muall, double deffall, char *namefile)
{
/**
 * prints results to file outside of working directory
 */
 
  char fnamemu[60];

  if(F >= 0){ 
  	sprintf(fnamemu, "../muoverfpos_R_%.2lf_setnumb_%d.dat", R_CONF, simparams.setnumb);
  }
  else{
  	sprintf(fnamemu, "../muoverfneg_R_%.2lf_setnumb_%d.dat", R_CONF, simparams.setnumb);
  }

  FILE *outmu;
  outmu=fopen(fnamemu, "a");
  fprintf(outmu, "%.3lf\t %.6lf\t %.6lf\t %s\n", F, muall/numtasks, deffall/numtasks, namefile);                       
  fclose (outmu); 

}

int histogramm_mpi_reduce(int m, 
                          int n, 
                          double backshift, 
                          double length, 
                          double bin, 
                          double **positions, 
                          char *fname, 
                          int taskid)
{
/**
 * function that stores 1 dimensional histogram in file named fname.
 * during a spatial sweep through all length/bin slices from position 
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
		for(k = 0; k < n; k++){
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
		  fprintf(outp,"%f\t %f\t %d\t %f\n", i*bin - backshift, (i+1)*bin - backshift, counterall, counterall/(simparams.N*bin));
		  fclose(outp);      
	  }
  }

  return countercheck;

}

int histogramm2d_mpi_reduce(int m, 
                            int n,
                            double bin2d, 
                            double **positionsx, 
                            double **positionsy,  
                            char *fname, 
                            int taskid)
{
/**
 * function that stores 2 dimensional histogram in file named fname.
 * during a spatial scan over all nx*ny rectangular fields, all particle
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
			for(j = 0; j < n; j++){
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
			  fprintf (outp, "%f\t%f\t%f\t%f\t\t%d\t\t%f\n", hx*bin2d, (hx+1)*bin2d, hy*bin2d - MAX_HALF_WIDTH, (hy+1)*bin2d - MAX_HALF_WIDTH, twodcounterall, twodcounterall/(simparams.N*bin2d*bin2d));
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
                       double initwidth, 
                       gsl_rng *r)
{
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
	  for(kset = 0; kset < simparams.setnumb; kset++){
		  do{
			  positionx[j][kset] = gsl_rng_uniform(r)*L;
			  positiony[j][kset] = (2*gsl_rng_uniform(r) - 1)*initwidth;
			    
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
          for(kset = 0; kset < simparams.setnumb; kset++){
	          fintx = 0;
		  finty = 0;
                  for(ktest = 0; ktest < simparams.setnumb; ktest++){
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

	if((xcheck != simparams.N) || (ycheck != simparams.N) || (twodcheck != simparams.N)){
	   printf("Error in Histogrammcounter!\n");
	}
}

double reset_pos_time(int setn_per_task, long int **posshift, long int **negshift) 
{
/**
 * Function to reset the sitcoeff.mulated system time t and the positions
 * to originate values in order to skip transient effects from
 * start of sitcoeff.mulation.
 */
	int i, j;
	double t = 0;
	for(i = 0; i < setn_per_task; i++){
	  for(j = 0; j < simparams.setnumb; j++){
		 posshift[i][j] = 0; 
		 negshift[i][j] = 0;

	  } 
	}
	return t;
}

void print_results_over_time(char *fname, 
                             char *fnamemom, 
                             double t, 
                             int abb, 
                             int abbdeff)
{
       /**
	* function to plot online results of sitcoeff.mulation over time 
	*/
	FILE *outp;
	outp=fopen(fname ,"a");
	fprintf(outp, "%.6f\t %.4Lf\t %.5lf\t %.4Lf\t %.5lf\t %d\t %d\n", t, tcoeff.meanx, tcoeff.mu, tcoeff.meanxsqu, tcoeff.deff, abb, abbdeff);
	fclose(outp);

	FILE *outpmom;
	outpmom=fopen(fnamemom ,"a");
	fprintf(outpmom, "%.6f\t %.6lf\t %.5Lf\n", t, tcoeff.meanspd, tcoeff.mthree);
	fclose(outpmom);
}

void adapt_posshifts(int shiftind, int i, int j, long int **posshift, long int **negshift)
{
/**
*function to update shifts which are monitored to calculate 
*absolute position in x-direction from sitcoeff.mulation with cyclic
*boundary conditions
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
 *function for the update of the counter that is used to monitore the
 *equilibration of the system. A transport quantity tran_quant is
 *compared with a previous value. If the difference is  below the    
 *demanded accurarcy, the value of the counter is increased.
 *If the difference is larger than twice the accurarcy the counter
 *is decreased
*/ 
	if(fabs(tran_quant - tran_quanto) <= accurarcy){
	  equcounter++;
	}

	if(fabs(tran_quant - tran_quanto) >= 2*accurarcy){
	  equcounter--;
	  if(equcounter < 0) equcounter = 0;
	}
	return equcounter;
}

int main (int argc, char **argv){
/** main function of Brownian motion sitcoeff.mulation
 *
 */

  /*
   * initialize parameters for time measurement and measure time
   */ 
  clock_t prgstart; 
  prgstart = clock(); 
  
  double u,v, x, y, yo, xo, yue, distx, disty, dist, f_cut,  deffall, meanspdall, muall;
  double muabbo, deffabbo;
  double t, dt, xtest, ytest;
  unsigned int j, xcountercheck, ycountercheck, twodcountercheck;
  long long i, testab, plotpoints, eq_stepnumbs;
  long double  meanxall, meanxsquall, msdall, mthreeall;
  long double  sqrt_flucts, f_dt;
  int abb, abbdeff, ktest, kset, shiftind; 
  int numtasks = 1;
  int taskid = MASTER;
  int setn_per_task;
//  int mpi_ind;
  char *namefile;
  char *conprfx;
  char *intprfx;
  struct par_specs *specs_args; 
  bool PosValid;
//  bool MPI_ON;

  /*
   * width of bins used for spatial discretization for position histograms   
   */
  double binx;
  double biny;
  double bin2d; 
  
  init_simparams(); 
  printf("\nN_test: %d\n", simparams.N);

  binx = 0.02*L;
  biny = 2.0*binx*MAX_HALF_WIDTH/L;
  bin2d = 0.05; 

  /*
   * read number of interacting particles per set from command line arguments
   */ 
  printf("\nargc: %d\n\n", argc);
  if(argc > 3){
	F = atof(argv[1]);
	simparams.setnumb = atof(argv[2]);
/*	mpi_ind = atof(argv[3]);
        if(mpi_ind = 1){
		MPI_ON = true;
        }
	else{
		MPI_ON = false;
	}*/
  }
  else{
 	return -1;
  }

  /*
   * get name of confinement and type of intra particle interaction
   * and provide them to function that creates working directory
   */
  conprfx = prfx_conf();
  intprfx = prfx_int();
  namefile = makedirectory(F, simparams.setnumb, conprfx, intprfx);


  /*
   * initialize arrays where x- and y-coordinates of particles are stored in
   */ 
  printf("\ninitialize arrays for positions and interactions\n"); 
  double **positionx;
  double **positiony;
  double **fintxarray;
  double **fintyarray;
  double **xstart;

  positionx = calloc(simparams.N/simparams.setnumb, sizeof(double));
  positiony = calloc(simparams.N/simparams.setnumb, sizeof(double));
  fintxarray = calloc(simparams.N/simparams.setnumb, sizeof(double));
  fintyarray = calloc(simparams.N/simparams.setnumb, sizeof(double));
  xstart = calloc(simparams.N/simparams.setnumb, sizeof(double));
  for(i = 0; i < simparams.N/simparams.setnumb; i++){
  	positionx[i] = calloc(simparams.setnumb, sizeof(double));
  	positiony[i] = calloc(simparams.setnumb, sizeof(double));
 	fintxarray[i] = calloc(simparams.setnumb, sizeof(double));
  	fintyarray[i] = calloc(simparams.setnumb, sizeof(double));
  	xstart[i] = calloc(simparams.setnumb, sizeof(double));
  }
  printf("arrays for positions and interactions are initialized!\n"); 


  long int **negshift;
  long int **posshift;
  negshift = calloc(simparams.N/simparams.setnumb, sizeof(long int));
  posshift = calloc(simparams.N/simparams.setnumb, sizeof(long int));
  for(i = 0; i < simparams.N/simparams.setnumb; i++){
  	negshift[i] = calloc(simparams.setnumb, sizeof(long int));
  	posshift[i] = calloc(simparams.setnumb, sizeof(long int));
  }
  

  char fname [60];
  char fnamemom [60];
  char fnamex [60];
  char fnamey [60];
  char fname2d [60];
  char fnamespecs [60];
  
#  ifdef MPI_ON
	  MPI_Init (&argc,&argv);
	  MPI_Comm_size (MPI_COMM_WORLD, &numtasks);
	  MPI_Comm_rank (MPI_COMM_WORLD, &taskid); 
# endif
 
  /**
   * calculate number of samples of interacting particles per task
   */ 
  setn_per_task = (int) simparams.N/(simparams.setnumb*numtasks);

  dt = time_step(B, R_INT, F); 
  sqrt_flucts = sqrt(2*BOTTRAD*dt);
  f_dt = F*dt;

  if(fabs(F) <= 0.1) n = simlong*n;

  if(F <= -2*pow(10,4)) n = 0.5*simlong*n;


  plotpoints = (long long) n/50; 
  testab = plotpoints*1;
  eq_stepnumbs = 30*testab;
  
  if(testab % plotpoints != 0){ 
	  printf("testab modulo plotpoints \n");
	  return -1;
  }  
  
  
  if(taskid == MASTER){
	    
	  sprintf(fname, "muovert_F_%.3lf.dat", F);
	  FILE *outp;
	  outp = fopen(fname ,"w");
	  fprintf(outp, "#time\t tcoeff.meanx\t tcoeff.mu\t Meansqdist\t  tcoeff.deff  \t abb\t abbdeff\n");
	  fclose(outp);
	  
	  sprintf(fnamemom, "momsovert_F_%.3lf.dat", F);
	  FILE *outpmom;
	  outpmom=fopen(fnamemom ,"w");
	  fprintf(outpmom, "#time\t tcoeff.meanspd\t  <x^3>-3<x^2><x>+2<x>^3\n");
	  fclose(outpmom);

	  if(simparams.N % (numtasks*simparams.setnumb) != 0){ 
		  printf("N modulo numtasks*setnumb\n");
		  return -1;
	  }

	  copy_main();
	  copycode_par(); 
	  copycode_conf();
	  copycode_int();
	  specs_args = par(n, dt, numtasks, testab, plotpoints); 
	  
          sprintf(fnamespecs, "muovert_specs.dat");
	  specs_basic(specs_args, fnamespecs);
	  specs_conf(binx, biny, bin2d);
	  f_cut = intforce(INT_CUTOFF,INT_CUTOFF);
	  specs_int(f_cut);

  	  if(numtasks > 1){
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
  
  init_particle_pos(setn_per_task, positionx, positiony, xstart, initwidth, r);

  printf("positions initialized!\n"); 

  /* Initialize inter-particle forces */
  init_particle_int(setn_per_task, positionx, positiony, fintxarray, fintyarray);

  printf("\ntask ID:\t %d -> particle positions and forces fixed\n", taskid);
  //printf("\nepsilon:\t %lf\n", EPS_L);
   
 
  /* Perform sitcoeff.mulation steps until equilibration is reached */
  t = 0;	
  i = 1;      
  abb = 0;	
  abbdeff = 0;	
  muabbo = 0;	
  deffabbo = 0;	
  do{	
	
	  t += dt;
	  i++;
	  /* Loop over trajectories */
	  for (j = 0; j < setn_per_task; j++){	
		  for(kset = 0; kset < simparams.setnumb; kset++){


			  xo = positionx[j][kset];
			  yo = positiony[j][kset];
			  
			  fintx = fintxarray[j][kset];
			  finty = fintyarray[j][kset];
		  
		 	  /* 
 			   * Perform sitcoeff.mulation steps until valid step where no particles overlap
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
					  /*sitcoeff.mulate particle-particle interaction*/ 
 					  if (simparams.setnumb > 1){ 
						    ktest = 0;
						    fintx = 0;
						    finty = 0;	
						    while((PosValid == true) && (ktest < simparams.setnumb)){
							 
							  xtest = positionx[j][ktest]; 
							  ytest = positiony[j][ktest];
	 
							  if(ktest != kset){
								  distx = x - xtest;
								  disty = y - ytest;

								  /* 
                                                                   * search relevant distance according
 								   * to minitcoeff.mum image conversion 
 								   */
								  if (abs(distx) > 0.5*L){ 
									distx = distx - L*(distx/abs(distx));
								  }
								  dist = sqrt(distx*distx + disty*disty);
								  if(dist <= 2*R_INT) PosValid = false;
								  
								  if((dist <= INT_CUTOFF) && (PosValid == true)){
									  fintxpair = intforce(distx, dist);
									  fintypair = intforce(disty, dist);
									  fintx += fintxpair;
									  finty += fintypair;

								  }
							  }

						    ktest++;
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
	
       	
	  if((i > n-testab) && (i <= n-testab+1)){ 
		  muabbo = tcoeff.mu;
		  deffabbo = tcoeff.deff;
	  }
	
          /*
           * Test progress of equilibration and plot results at certain 
           * sitcoeff.mulation steps i
           */ 	
	  if(((i > n) && (i%testab == 0)) || (i%plotpoints == 0)){ 

          /*call function for calculation of transport coefficients 
	   * such as mobility or mean-squared displacement */ 
          calc_transpcoeffs(numtasks,
	       	            setn_per_task, 
		            t, F, 
		            posshift,
		            negshift,
		            positionx,
		            xstart);
	  
			 
                  /*
                   * Update of the equilibration counter that are used to monitore the
                   * equilibration of the mobility and diffusivity
                   */  
		  if((i > n) && (i % testab == 0)){ 
 			  abb = update_equcounter(tcoeff.mu, muabbo, accur, abb);

			  if(tcoeff.deff > 1.0){  
				  
 			  	  abbdeff = update_equcounter(tcoeff.deff, deffabbo, deffaccur*deffabbo, abbdeff);
			  }

			  else{ 
 			  	  abbdeff = update_equcounter(tcoeff.deff, deffabbo, deffaccur, abbdeff);
			  }		
			  muabbo = tcoeff.mu;    
			  deffabbo = tcoeff.deff;
                  }
 		  /*
                   * Plot results to check progress of equilibration 
                   */
		  if((i%plotpoints == 0) && (taskid == MASTER)){
			  print_results_over_time(fname, 
						  fnamemom, 
						  t, 
						  abb, 
						  abbdeff);

		  }
	  }

	  /*
           *  Reset position and time information to truncate transient effects from small times 
           */
	  if(i == eq_stepnumbs){
		  t = reset_pos_time(setn_per_task, posshift, negshift); 
	  } 
 
  /* 
   * Closes while loop over sitcoeff.mulation steps if criteria for equilibration are fulfilled
   */
  }while((abb < numbtest) || (abbdeff < numbtest));
 
  /*Merge quantities that were calculated in separated MPI threads*/ 
#  ifdef MPI_ON
	  MPI_Reduce(&tcoeff.meanx, &meanxall, 1, MPI_LONG_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
	  MPI_Reduce(&tcoeff.meanxsqu, &meanxsquall, 1, MPI_LONG_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
	  MPI_Reduce(&tcoeff.msd, &msdall, 1, MPI_LONG_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
	  MPI_Reduce(&tcoeff.mthree, &mthreeall, 1, MPI_LONG_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
	  MPI_Reduce(&tcoeff.meanspd, &meanspdall, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
	  MPI_Reduce(&tcoeff.mu, &muall, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
	  MPI_Reduce(&tcoeff.deff, &deffall, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
#  else
        meanxall = tcoeff.meanx;
        meanxsquall = tcoeff.meanxsqu;
        msdall = tcoeff.msd;
        mthreeall = tcoeff.mthree;
        meanspdall = tcoeff.meanspd;
        muall = tcoeff.mu;
        deffall = tcoeff.deff;
#  endif
  
  print_runtime_threads(prgstart, numtasks, taskid);
  
  sprintf(fnamex, "meanx_Histogram_F_%.3lf.dat", F);
  xcountercheck = histogramm_mpi_reduce(setn_per_task, 
                                        simparams.setnumb, 
					0, 
					L, 
                                        binx, positionx, 
                                        fnamex, taskid);

  sprintf(fnamey, "meany_Histogram_F_%.3lf.dat", F);
  ycountercheck = histogramm_mpi_reduce(setn_per_task, 
                                        simparams.setnumb, 
                                        MAX_HALF_WIDTH, 
                                        2*MAX_HALF_WIDTH, 
                                        biny, positiony, 
                                        fnamey, 
                                        taskid);

  sprintf(fname2d, "meanpos_Histogram2d_F_%.3lf.dat", F);
  twodcountercheck =  histogramm2d_mpi_reduce(setn_per_task, 
                                              simparams.setnumb, bin2d, 
                                              positionx, positiony, 
                                              fname2d, taskid);


  if(taskid == MASTER){
 	   print_positions(setn_per_task, simparams.setnumb, positionx, positiony);
           print_hist_countercheck(xcountercheck, ycountercheck, twodcountercheck, fnamespecs);
           print_resallthreads(msdall, meanspdall, muall, deffall, meanxall, meanxsquall, mthreeall, numtasks, &fname[0], &fnamemom[0]);
           print_muoverf(F, numtasks, muall, deffall, namefile);
            
           delerrorfiles(F,simparams.setnumb);           
            
	   print_runtime(prgstart, numtasks);
           
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

  return 0;
}
