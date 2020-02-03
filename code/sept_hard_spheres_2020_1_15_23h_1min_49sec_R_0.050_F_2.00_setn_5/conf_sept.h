#define B 0.1
#define AMP 1.0
#define MAX_HALF_WIDTH (AMP+B)
#define R_CONF_SQ R_CONF*R_CONF

extern double yuef_ext(double x, double y);
extern double yuef_sept(double x, double y);
extern void specs_conf(double binx, double biny, double bin2d);
extern void copycode_conf();
extern char *prfx_conf(); 
