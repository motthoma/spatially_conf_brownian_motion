#include "print_routines.h"
#include "par_sim.h"
#include "results_transport.h"

#include <stdio.h>

void PRINT_results_over_time(double t, 
                             int abb, 
                             int abbdeff)
{
       /**
	* function to print online results of simulation over time 
	*/
/*	printf("start_to_print_result\n"); */
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


void PRINT_copycode(){

   char copycode[200];

  sprintf(copycode, "cp ../print_routines.* ./");
  system(copycode);

}
