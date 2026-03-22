# include "../src/conf_cos.h"

/**
 * Provide wrapper functions to expose yuef_conf
 * e.g. for plotting the confinement profile in python. 
 */
double CONF_yuef_cos(double x, double y) {
    return yuef_cos(x, y);
}

/**
 * Function provides boundary of confinement with cosine profile.
 *
 */
double CONF_yu_cos(double x){
	return (BOTTLENECK_WIDTH + AMP + AMP*sin(K_COS*x));
}

/**
 * Function provides effective boundary of confinement with cosine profile
 * for a particle with radius R_CONF.
 * It computes the maximal allowed y-position for the particle center at position x.
 */
double CONF_yu_eff_cos(double x){
    double xcirc;
    double xt;
    double yu;
    double y_bound;
    double y_min = 1e10;

    if(R_CONF == 0){
        return (BOTTLENECK_WIDTH + AMP + AMP*sin(K_COS*x));
    }

    /* Scan the particle's surface to find the minimum distance to the boundary */
    for(int i = 0; i <= CHECKN; i++){
        xcirc = -R_CONF + 2.0*R_CONF*i/CHECKN;
        xt = x + xcirc;
        yu = BOTTLENECK_WIDTH + AMP + AMP*sin(K_COS*xt);
        y_bound = yu - sqrt(R_CONF*R_CONF - xcirc*xcirc);
        if(y_bound < y_min){
            y_min = y_bound;
        }
    }
    return y_min;
}
