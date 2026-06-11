// sim_config.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h> // Required for string functions
#include <math.h>
#include <stdbool.h>
#include "sim_config.h"
#include "code_handling.h"
#include "comp_gen_header.h"

extern T_SimParams SimParams; // Declare SimParams as external

/**
 * Reads simulation parameters from a configuration file.
 * @param filepath The path to the configuration file.
 * @return true if parameters were read successfully, false otherwise.
 */
bool SIMCONFIG_read_params(const char *filepath) {
    FILE *file = fopen(filepath, "r");
    if (file == NULL) {
        fprintf(stderr, "Error: Could not open config file '%s'\n", filepath);
        return false;
    }

    char line[256]; // Buffer to hold each line
    while (fgets(line, sizeof(line), file) != NULL) {
        // Skip comments and empty lines
        if (line[0] == '#' || strlen(line) <= 1) {
            continue;
        }

        char *key = strtok(line, "=");
        char *value = strtok(NULL, "\n"); // Read until newline

        // Trim whitespace
        if (key) {
            // Trim leading whitespace
            while (*key == ' ' || *key == '\t') {
                key++;
            }
            // Trim trailing whitespace
            char *end = key + strlen(key) - 1;
            while (end > key && (*end == ' ' || *end == '\t')) {
                end--;
            }
            *(end + 1) = '\0';
        }

        if (value) {
            // Trim leading whitespace
            while (*value == ' ' || *value == '\t') {
                value++;
            }
            // Trim trailing whitespace
            char *end = value + strlen(value) - 1;
            while (end > value && (*end == ' ' || *end == '\t')) {
                end--;
            }
            *(end + 1) = '\0';
        }

        if (key != NULL && value != NULL) {
            if (strcmp(key, "F") == 0) {
                SimParams.F = atof(value);
            } else if (strcmp(key, "parts_per_set") == 0) {
                SimParams.parts_per_set = atoi(value);
            } else if (strcmp(key, "N") == 0) {
                SimParams.N = atoi(value);
            } else if (strcmp(key, "patience") == 0) {
                SimParams.patience = atoi(value);
            } else if (strcmp(key, "stepnumb") == 0) {
                SimParams.stepnumb = atoi(value);
            } else if (strcmp(key, "accur") == 0) {
                SimParams.accur = atof(value);
            } else if (strcmp(key, "initwidth") == 0) {
                SimParams.initwidth = atof(value);
            } else if (strcmp(key, "init_max_xpos") == 0) {
                SimParams.init_max_xpos = atof(value);
            }
            // Add more parameters here as needed
        }
    }

    fclose(file);
    return true;
}

T_SimParams SimParams;

/**
 * Initialize simulation parameters with default values.
 */
void SIMCONFIG_init(int tasks, const char *config_filepath) {
    /* Set basic parameters that might be needed before reading config, or are not configurable */
    SimParams.N = 100; // Default N, can be overridden by config
    SimParams.numtasks = tasks; // Set numtasks from argument

    /* Set internal defaults for all parameters. These can be overridden by the config file. */
    SimParams.patience = 20;
    SimParams.stepnumb = 50000;       // 5*1e4
    SimParams.accur = 1e-2;
    SimParams.initwidth = 1.0;
    SimParams.init_max_xpos = 0.7;
    SimParams.F = 0.0; // Default force
    SimParams.parts_per_set = 1; // Default parts per set

    // Attempt to read parameters from the config file.
    // Values read from the file will override the defaults set above.
    if (!SIMCONFIG_read_params(config_filepath)) {
        fprintf(stderr, "Warning: Could not read %s. Using internal defaults for F and parts_per_set.\n", config_filepath);
        // No need to set defaults for F and parts_per_set again here, they are already set above
    }
    
    /* Calculate derived parameters based on the current SimParams values (defaults or read from file) */
    SimParams.print_interval_steps = SimParams.stepnumb / 100;
    SimParams.testab = SimParams.print_interval_steps;
    SimParams.reset_stepnumb = 15 * SimParams.testab;
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
    int min_steps = SimParams.stepnumb + SimParams.testab * SimParams.patience;
    double eq_time = SimParams.reset_stepnumb * SimParams.time_step;
    double total_time = (SimParams.stepnumb - SimParams.reset_stepnumb 
                        + SimParams.testab * SimParams.patience) * SimParams.time_step;
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
            "Accuracy of mobility: %lf\n",
            SimParams.N,
            SimParams.parts_per_set,
            SimParams.accur);

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
            SimParams.patience,
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
    CODEHAND_copy_file_to_dest("sim_params.conf");
}
