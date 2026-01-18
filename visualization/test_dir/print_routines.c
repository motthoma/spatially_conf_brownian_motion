#include "print_routines.h"

/* Define the global variable exactly once */
struct PrintResults Print = {"", "", false};

#include "sim_config.h"
#include "code_handling.h"
#include "results_transport.h"

#include <stdio.h>
#include <time.h>

void PRINT_muoverf(double muall, double deffall)
{
/**
 * prints results to file outside of working directory
 */
  char fnamemu[60];
  snprintf(fnamemu, 
           sizeof fnamemu,
           "%s/muoverfpos_R_%.2lf_parts_per_set_%d.dat",
           RUNS_DIR,
           R_CONF,
           SimParams.parts_per_set);
  FILE *outmu;
  outmu=fopen(fnamemu, "a");
  fprintf(outmu, "%.3lf\t %.6lf\t %.6lf\t %s\n",
                  SimParams.F,
                  muall/SimParams.numtasks,
                  deffall/SimParams.numtasks,
                  DestPaths.destdir_name);
  fclose (outmu); 
}

void PRINT_set_print_flag(int time_step)
{
/*
 * Test progress of equilibration and plot results at certain
 * simulation steps i
 */

  Print.PrintRes = false;
  if(time_step % SimParams.print_interval_steps == 0){
      Print.PrintRes = true;
  }
}

void PRINT_header_for_results_over_time()

{
    /**
    * function to print online results of simulation over time 
    */
    snprintf(Print.fname,
             sizeof Print.fname,
             "%s/muovert_F_%.3lf.dat",
             DestPaths.fullpath,
             SimParams.F);

    FILE *outp;
    outp = fopen(Print.fname ,"w");
    fprintf(outp, "#time\t meanx: <x>\t meanspeed: <v>\t mu\t abb\n");
    fclose(outp);

    printf("#time\t meanx: <x>\t meanspeed: <v>\t mu\t abb\n");

    snprintf(Print.fnamemom,
             sizeof Print.fnamemom,
             "%s/momsovert_F_%.3lf.dat",
             DestPaths.fullpath,
             SimParams.F);
    FILE *outpmom;
    outpmom=fopen(Print.fnamemom ,"w");
    fprintf(outpmom, "#time\t meanx: <x>\t meanxsquare: <x^2>\t Meansqdist: <x^2> - <x>^2\t  deff: (<x^2> - <x>^2)/(2t)\t third cumulant: <x^3>-3<x^2><x>+2<x>^3\n");
    fclose(outpmom);
}

void PRINT_results_over_time(double t, 
                             int abb) 
{
/*
* function to print online results of simulation over time 
*/

    if(Print.PrintRes == true){
        FILE *outp;
        outp=fopen(Print.fname ,"a");
        fprintf(outp, "%.6f\t %.4Lf\t %.5lf\t %.4lf\t  %d\n",
                      t,
                      tcoeff.meanx,
                      tcoeff.meanspeed,
                      tcoeff.mu,
                      abb);
        fclose(outp);
        printf("%.6f\t %.4Lf\t %.5lf\t %.4lf\t  %d\n",
               t,
               tcoeff.meanx,
               tcoeff.meanspeed,
               tcoeff.mu,
               abb);

        FILE *outpmom;
        outpmom=fopen(Print.fnamemom ,"a");
        fprintf(outpmom, "%.6f\t %.6Lf\t %.5Lf\t %.5Lf\t %.5lf\t %.5Lf\n",
                         t,
                         tcoeff.meanx,
                         tcoeff.meanxsqu,
                         tcoeff.msd,
                         tcoeff.deff,
                         tcoeff.thirdcum);
        fclose(outpmom);
    }
}

void PRINT_positions(double **posx, double **posy){
/**
 * prints particle positions to file
 */
  int i;
  int j;
  FILE *outpos;

  char positions[256];
  snprintf(positions, 256, "%s/positions.dat", DestPaths.fullpath); 
  outpos=fopen(positions, "a");
  fprintf(outpos, "#xpositions\t ypositions\n");

  for(i = 0; i < SimParams.setn_per_task; i++){
	for(j = 0; j < SimParams.parts_per_set; j++){
		fprintf(outpos, "%.5lf\t %.5lf\n", posx[i][j], posy[i][j]);
	}

  }
  fclose(outpos);
}

void PRINT_runtime(clock_t start)
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
  outpspecs=fopen(DestPaths.fname_simparams, "a");
  fprintf(outpspecs, "\nComputing time: %d days %dhours %dmin %dsec\n",
          timediff/(3600*24),
          (timediff/3600)%24, 
          (timediff/60)%60,
          timediff%60);

  fprintf(outpspecs, "Total computing time of all threads: %d days %dhours %dmin %dsec\n\n",
          timediff_all/(3600*24), 
          (timediff_all/3600)%24, 
          (timediff_all/60)%60, 
          timediff_all%60);
  fclose(outpspecs);

}


void PRINT_runtime_threads(clock_t start,
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
    char taskres [256];
    snprintf(taskres, 256, "%s/taskres.dat", DestPaths.fullpath); 
	outptasks=fopen(taskres, "a");
	fprintf(outptasks, "\nTask %d:\nmeanx = %.8Lf\t meanxsqu = %.8Lf\t meanspeed = %.8lf\t mu = %.8lf\t deff = %.8lf\n",
            taskid, 
            tcoeff.meanx, 
            tcoeff.meanxsqu,
            tcoeff.meanspeed, 
            tcoeff.mu, 
            tcoeff.deff);

	fprintf(outptasks, "Time of simulation for task %d: %d days %dhours %dmin %dsec\n",
            taskid, 
            timediff/(3600*24), 
            (timediff/3600)%24, 
            (timediff/60)%60, 
            timediff%60);
	fclose(outptasks);
  }
}

void PRINT_resallthreads(long double msdall, 
                         double meanspeedall, 
                         double muall, 
                         double deffall, 
                         long double meanxall, 
                         long double meanxsquall, 
                         long double thirdcumall)
{
  /**
  * prints results of simulation to file at the end of the simulation
  */
  FILE *outp;
  outp=fopen(Print.fname ,"a");
  fprintf(outp, "\n\nAverage of all Threads:\n\nmeanx = %.5Lf\t meanspeed = %.5lf\t mu = %.5lf\n\n", 
          meanxall/SimParams.numtasks, 
          meanspeedall/SimParams.numtasks, 
          muall/SimParams.numtasks);
  fclose(outp);
           
  FILE *outpmom;
  outpmom=fopen(Print.fnamemom ,"a");
  fprintf(outpmom, "\n\nAverage of all Threads:\n\nmeanx = %.5Lf\t meanxsqu = %.5Lf\t msdall = %.5Lf\t deff = %.5lf\t thirdcum = %.5LF\n\n", 
          meanxall/SimParams.numtasks, 
          meanxsquall/SimParams.numtasks, 
          msdall/SimParams.numtasks,
          deffall/SimParams.numtasks,  
          thirdcumall/SimParams.numtasks);
  fclose(outpmom);
}



void PRINT_copycode(){
    CODEHAND_copy_file_to_dest("print_routines.c");
    CODEHAND_copy_file_to_dest("print_routines.h");
}
