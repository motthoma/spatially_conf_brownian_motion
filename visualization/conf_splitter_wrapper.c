# include "../src/conf_splitter.h"

/**
 * Provide wrapper functions to expose yuef_splitter
 * e.g. for plotting the confinement profile in python. 
 */
double CONF_yuef_wrapper(double x, double y) {
    return yuef_splitter(x, y);
}

/**
 * Function provides boundary of confinement with saw-tooth profile
 * used for the entropic splitter (see Motz et al. J. Chem. Phys. 2014).
 * Function is currently only used for plotting the confinement 
 * profile in python, but can be used for other purposes as well.
 *
 * Boundary can only be represented in one region: At the
 * bottleneck, yu has a jump from the bottleneck to the widest part.
 * From the widest point, a straight line with slope M gives a
 * decending boundary.
 */
static inline double yu_splitter(double x, double y){
    if(x == 0){
        return BOTTLENECK_WIDTH;
    }
    else{
        return (BOTTLENECK_WIDTH + SLOPE*(L_CONF - x));
    }
}
