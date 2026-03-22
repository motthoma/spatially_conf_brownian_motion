// sim_config.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "sim_config.h"
#include "code_handling.h"
#include "comp_gen_header.h"

T_SimParams SimParams;

/**
 * Initialize simulation parameters with default values.
 */
void SIMCONFIG_init(int tasks) {
    /* user config parameters*/
    SimParams.N = 1000;
    SimParams.numbtest = 20;
    SimParams.stepnumb = 50000;       // 5*1e4
    SimParams.simlong = 20;
    SimParams.accur = 1e-2;
    SimParams.deffaccur = 1e7;
    SimParams.initwidth = 1.0;
    SimParams.init_max_xpos = 0.7;

    SimParams.print_interval_steps = SimParams.stepnumb / 100;
    SimParams.testab = SimParams.print_interval_steps;
    SimParams.reset_stepnumb = 15 * SimParams.testab;
    SimParams.numtasks = tasks;
    SimParams.n_interact_sets = (int) SimParams.N/SimParams.parts_per_set;
    /* number of particle ensembles simulated in one task.
     * this number has to be at least one, to prevent 
     * necessary communication between the threads during 
     * simulation time.
     */
    SimParams.setn_per_task = (int) SimParams.N/(SimParams.parts_per_set*SimParams.numtasks);
    /* calculate number of sets of interacting particles */ 
    SimParams.time_step = SIMCONFIG_time_step(BOTTLENECK_WIDTH, R_INT); 
}

/**
 * Check consistency of simulation parameters.
 * Returns true if consistent, false otherwise.
 */
bool SIMCONFIG_check_consistency() {
    bool consistent = true;

    if (SimParams.testab % SimParams.print_interval_steps != 0) {
        printf("Error: testab modulo print_interval_steps != 0\n");
        consistent = false;
    }

    if (SimParams.N % (SimParams.numtasks * SimParams.parts_per_set) != 0) {
        printf("Error: N modulo (numtasks * parts_per_set) != 0\n");
        consistent = false;
    }

    if (!consistent) {
        printf("Chosen simulation parameters are inconsistent!\n");
    }

    return consistent;
}

/**
 * Calculate simulation time step based on confinement and particle scales.
 */
double SIMCONFIG_time_step(double lscale_conf, double lscale_part) {
    double tstep_scale = (SimParams.parts_per_set == 1) ? lscale_conf : lscale_part;
    double tstep;

    if (fabs(SimParams.F) <= 1 / lscale_conf) {
        tstep = tstep_scale * tstep_scale * TSTEP_BASE;
    } else {
        tstep = (tstep_scale / fabs(SimParams.F)) * TSTEP_BASE;
    }

    return tstep * L_CONF;
}

/**
 * Write simulation parameters to a specs file for documentation.
 */
void SIMCONFIG_write_specs() {
    int min_steps = SimParams.stepnumb + SimParams.testab * SimParams.numbtest;
    double eq_time = SimParams.reset_stepnumb * SimParams.time_step;
    double total_time = (SimParams.stepnumb - SimParams.reset_stepnumb 
                        + SimParams.testab * SimParams.numbtest) * SimParams.time_step;
    double check_time = SimParams.testab * SimParams.time_step;
    double readout_time = SimParams.print_interval_steps * SimParams.time_step;
  
    char fname_simparams [60]; 
    snprintf(fname_simparams,
             sizeof fname_simparams,
             "parameters_simulation_overall.dat");
    DestPaths.fname_simparams = malloc(256);    
    snprintf(DestPaths.fname_simparams,
             800,
             "%s/%s",
             DestPaths.fullpath, fname_simparams);

    FILE *out_file = fopen(DestPaths.fname_simparams, "w");
    if (!out_file) {
        return;
    }

    fprintf(out_file, "\nSimulation Parameters:\n\n");

    fprintf(out_file,
            "# of particles: %d\n"
            "# of particles per set: %d\n\n"
            "Accuracy of mobility: %lf\n"
            "Accuracy of Deff: %lf\n\n",
            SimParams.N,
            SimParams.parts_per_set,
            SimParams.accur,
            SimParams.deffaccur);

    fprintf(out_file,
            "Number of steps until equilibration checks start: %.2e\n"
            "Time step size for Langevin propagation: %.2e\n"
            "Time until equilibration is reset: %.2lf\n"
            "Number of equilibration checks: %d\n"
            "Steps between two tests: %.2e\n"
            "Steps between readouts: %.2e\n"
            "Time between accuracy checks: %.2lf\n"
            "Time between readouts: %.2lf\n"
            "Minimum number of simulation steps: %.2e\n"
            "Minimum simulation time: %.3lf\n\n",
            (double)SimParams.stepnumb,
            SimParams.time_step,
            eq_time,
            SimParams.numbtest,
            (double)SimParams.testab,
            (double)SimParams.print_interval_steps,
            check_time,
            readout_time,
            (double)min_steps,
            total_time);

    fprintf(out_file,
            "# of parallelized tasks: %d\n"
            "Applied force F: %.2lf\n"
            "Width of initial particle distribution: %.2lf\n"
            "Maximum initial x-coordinate: %.2lf\n\n",
            SimParams.numtasks,
            SimParams.F,
            SimParams.initwidth,
            SimParams.init_max_xpos);

    fclose(out_file);
}

/**
 * Copy the sim_config module to the working directory.
 */
void SIMCONFIG_copy_code() {
    CODEHAND_copy_file_to_dest("sim_config.c");
    CODEHAND_copy_file_to_dest("sim_config.h");
}
