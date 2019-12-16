#define AMP (1.0/(2.0*M_PI))
#define MAX_HALF_WIDTH (2*AMP+B)

#define CHECKN 5
#define K_COS (2.0*M_PI/L)

#define XMAX (L/4.0-R_CONF)
#define SQRT_SHIFT R_CONF*sqrt(1+AMP*AMP*K_COS*K_COS)

extern double yuef_ext(double x, double y);
extern double yuef_cos(double x, double y);
extern void specs_conf(double binx, double biny, double bin2d);
extern void copycode_conf();
extern char *prfx_conf(); 
