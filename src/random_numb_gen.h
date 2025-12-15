#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


/**
 *********************************************************
 *
 * External types and variables
 *
 *********************************************************
 */

// Struct to hold core simulation values
typedef struct {
    gsl_rng *r;  // Pointer to GSL random number generator
} T_GSL_RNG;
extern T_GSL_RNG GSL_RNG;


/**
 * @brief Initialize the gsl_rng pointer inside SimCoreVals
 *
 * @param sim Pointer to SimCoreVals struct
 * @param taskid Integer used to vary the seed
 */
void SIM_init_rng(int taskid, unsigned long seed);

/**
 * @brief Get a normally distributed random number
 * @param sigma Standard deviation (σ)
 * @return Gaussian random number with mean 0 and stddev σ
 */
double SIM_get_gaussian(double mu, double sigma);

/**
 * @brief Get a uniform random number in [0,1)
 */
double SIM_get_uniform(void);
    
/**
 * @brief Free the global GSL random number generator.
 *
 * Frees the allocated RNG, if any, and sets the pointer to NULL.
 */
void SIM_free_rng(void); 
