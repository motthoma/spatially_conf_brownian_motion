#define R_INT (1.0*R_CONF)
#define INT_CUTOFF 0

extern double fintx, finty, fintxpair, fintypair;

static inline double intforce(double disti, double dist){

	return(0);

}
//extern int intforce(double disti, double dist);
extern void specs_int(double f_cut);
extern void copycode_int();
extern char *prfx_int(); 
