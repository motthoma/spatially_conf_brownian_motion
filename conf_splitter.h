#define M 0.9
#define MAX_HALF_WIDTH (M*L+B)

#define SQRT_SHIFT R_CONF*sqrt(1+M*M)
#define R_CONF_SQ  R_CONF*R_CONF
#define Lp (L-(R_CONF*M)/(sqrt(1+M*M)))


extern double yuef_ext(double x, double y);
extern double yuef_splitter(double x, double y);
extern void specs_conf(double binx, double biny, double bin2d);
extern void copycode_conf();
extern char *prfx_conf(); 
