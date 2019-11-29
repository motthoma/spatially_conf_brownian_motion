#define R_INT (R_CONF)
#define INT_CUTOFF (0.2*L)
#define EPS_L (1.0)
#define LJMIN (0.5*B)
#define LJMINPOW (pow(LJMIN,6))
#define LJPREFAC 12.0*EPS_L*LJMINPOW

extern double fintx, finty, fintxpair, fintypair;

extern double intforce(double disti, double dist);
extern double intforce_LJ(double disti, double dist);
extern void specs_int(double f_cut);
extern void copycode_int();
extern char *prfx_int(); 
