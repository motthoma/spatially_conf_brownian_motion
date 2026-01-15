#include <math.h>
#include "equilibration_manager.h"



void EQUIMAN_init(T_EquManager *EquManager, int stepnumb, int testab)
{
    EquManager->equ_counter = 0;
    EquManager->mu_old = 0;
    EquManager->test_start_step = stepnumb - testab;
}

void EQUIMAN_update_mu_old(T_EquManager *EquManager, int time_step, double mu)
{
/*
 *Initialize reference values of mobility and diffusion coefficients
 *for later judgement of equilibration process.
 */

    if((time_step > EquManager->test_start_step) && (time_step <= EquManager->test_start_step + 1)){ 
      EquManager->mu_old = mu;
    }
}

void EQUIMAN_set_test_flag(T_EquManager *EquManager, int time_step, int stepnumb, int testab)
{
/*
 * Test progress of equilibration and plot results at certain
 * simulation steps i
 */

  EquManager->TestRes = false;
  if((time_step > stepnumb) && (time_step % testab == 0)){
      EquManager->TestRes = true;
  }
}

void EQUIMAN_update_counter(T_EquManager *EquManager, double mu_current,
                               double accurarcy)
{
/*
 * Function for the update of the counter that is used to monitore the
 * equilibration of the system. A transport quantity tran_quant is
 * compared with a previous value. If the difference is  below the    
 * demanded accurarcy, the value of the counter is increased.
 * If the difference is larger than the accurarcy the counter
 * is decreased
 */ 
    if (EquManager->TestRes == true){
        if(fabs(mu_current - EquManager->mu_old) <= accurarcy){
          EquManager->equ_counter++;
        }
        else{
          EquManager->equ_counter--;
          if(EquManager->equ_counter < 0){
                EquManager->equ_counter = 0;
          }
        }
        EquManager->mu_old = mu_current;
    }
}
