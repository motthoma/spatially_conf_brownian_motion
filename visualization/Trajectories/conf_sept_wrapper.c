# include "../src/conf_sept.h"

/**
 * Provide wrapper functions to expose yuef_sept
 * e.g. for plotting the confinement profile in python. 
 */
double CONF_yuef_sept(double x, double y) {
    return yuef_sept(x, y);
}

/**
 * Function provides boundary of confinement with septated channel.
 */
double CONF_yu_sept(double x){
	return (BOTTLENECK_WIDTH + AMP - R_CONF);
}
