#include <stdio.h>
#include <time.h>

#include "comp_gen_header.h"
#include "code_handling.h"


char* CODEHAND_makedirectory(char *confprfx, char *intprfx){             
/**
 * creates directory, where code and data is transferred to
 * moves to this directory as working directory
 */
  time_t tnow;
  struct tm *tmnow;
  time(&tnow);
  tmnow = localtime(&tnow);

  char *e = malloc(100*sizeof(char));
  #ifdef LJMIN
	  sprintf(e, "%s_%s_%d_%d_%d_%dh_%dmin_%dsec_R_%.3lf_LJMIN_%.3lf_EPS_%.2lf_F_%.2lf_setn_%d", confprfx, 
												     intprfx, 
												     tmnow->tm_year + 1900, 
												     tmnow->tm_mon+1, 
												     tmnow->tm_mday, 
												     tmnow->tm_hour,
												     tmnow->tm_min,
												     tmnow->tm_sec, 
												     R_CONF,
												     LJMIN,
												     EPS_L, 
												     SimParams.F, 
												     SimParams.setnumb);
  #else
	 
	  sprintf(e, "%s_%s_%d_%d_%d_%dh_%dmin_%dsec_R_%.3lf_F_%.2lf_setn_%d", confprfx, 
									       intprfx, 
									       tmnow->tm_year + 1900, 
									       tmnow->tm_mon+1, 
									       tmnow->tm_mday, 
									       tmnow->tm_hour,
									       tmnow->tm_min,
									       tmnow->tm_sec, 
									       R_CONF, 
									       SimParams.F, 
									       SimParams.setnumb);
  #endif
  /**
   * 0700 is modus for determining access rights to dir
   */
  mkdir(e, 0700);
  chdir(e);
  return e;
}

void CODEHAND_copy_main(){
/**
 * copies main_brownianconf.c to directory created by 'makedirectory'
 */

  char copycode[200];

  sprintf(copycode, "cp ../main_brownconf.c ../masterinteract.py ../makefile ./");
  system(copycode);

}

void CODEHAND_delerrorfiles(){
/**
 * deletes error and log files used on albeniz
 */
  char delerrorfile[100];
  sprintf(delerrorfile, "rm ../F_%.2lf_sn_%.0d_clintparallel* ",SimParams.F, SimParams.setnumb);
  system(delerrorfile);

}

void CODEHAND_copycode(){

   char copycode[200];

  sprintf(copycode, "cp ../code_handling.* ./");
  system(copycode);

}


void CODEHAND_copy_comp_gen_header(){
 			char copycode[200];
 			sprintf(copycode, "cp ../comp_gen_header.h ./");
 			system(copycode);
}
