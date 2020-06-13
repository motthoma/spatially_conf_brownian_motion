#include <stdio.h>

#include "par_sim.h"
#include "comp_gen_header.h"
#include "results_transport.h"


T_TransportCoeffs tcoeff;


void TCOEF_init(){
	
        /*average of particle x-positions <x>*/
	tcoeff.meanx = 0; 
	/*average of particle speed in x-directions <v>*/
	tcoeff.meanspeed = 0.0;
	/*<x^2>*/
	tcoeff.meanxsqu = 0.0; 
	/*<x^3>*/
        tcoeff.meanxqub = 0.0;
        /*mean-squared displacement of x-position (2nd moment): <x^2> - <x>^2*/
        tcoeff.msd = 0.0;
	/*third cumulant of x-position: <x^3> - 3*<x>*<x^2> + 2*<x>^3*/
	tcoeff.thirdcum = 0.0;
        /*effective diffusion coeff.: (<x^2> - <x>^2)/(2*t*B/R)*/
	tcoeff.deff = 0.0; 
	/*non-linear mobility mu = <v>/F*/
	tcoeff.mu = 0.0;

}

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
