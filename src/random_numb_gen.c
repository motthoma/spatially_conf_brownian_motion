#include "random_numb_gen.h"
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <gsl/gsl_rng.h>

T_GSL_RNG GSL_RNG;


/**
 * @brief Initialize the global GSL random number generator
 * @param taskid Integer used to vary the seed
 */

void RNG_init_rng(int taskid, unsigned long seed) {
    if (GSL_RNG.r) return; // safety check

    // Allocate RNG
    GSL_RNG.r = gsl_rng_alloc(gsl_rng_mt19937);
    if (!GSL_RNG.r){ 
        printf("random number generator initialization failed!");
        return; // safety check
    }

    // Seed RNG with current time + taskid
    gsl_rng_set(GSL_RNG.r, seed);
}

/**
 * @brief Get a normally distributed random number
 * @param sigma Standard deviation (σ)
 * @return Gaussian random number with mean 0 and stddev σ
 */
double RNG_get_gaussian(double mu, double sigma) {
    assert(GSL_RNG.r && "RNG not initialized");
    return mu + gsl_ran_gaussian_ziggurat(GSL_RNG.r, sigma);
}

/**
 * @brief Get a uniform random number in [0,1)
 */
double RNG_get_uniform(void) {
   assert(GSL_RNG.r && "RNG not initialized");
   return  gsl_rng_uniform(GSL_RNG.r);
}

/**
 * @brief Free the global GSL random number generator.
 *
 * Frees the allocated RNG, if any, and sets the pointer to NULL.
 */
void RNG_free_rng(void) {
    if (GSL_RNG.r) {
        gsl_rng_free(GSL_RNG.r);
        GSL_RNG.r = NULL;
    }
}
