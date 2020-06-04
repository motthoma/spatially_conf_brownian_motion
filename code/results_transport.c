#include "results_transport.h"


T_TransportCoeffs tcoeff;


void TCOEF_init(){
	
        /*average of particle x-positions <x>*/
	tcoeff.meanx = 0; 
	/*average of particle speed in x-directions <v>*/
	tcoeff.meanspeed = 0.0;
	/*<x^2>*/
	tcoeff.meanxsqu = 0.0; 
	/*<x^3>*/
        tcoeff.meanxqub = 0.0;
        /*mean-squared displacement of x-position (2nd moment): <x^2> - <x>^2*/
        tcoeff.msd = 0.0;
	/*third cumulant of x-position: <x^3> - 3*<x>*<x^2> + 2*<x>^3*/
	tcoeff.thirdcum = 0.0;
        /*effective diffusion coeff.: (<x^2> - <x>^2)/(2*t*B/R)*/
	tcoeff.deff = 0.0; 
	/*non-linear mobility mu = <v>/F*/
	tcoeff.mu = 0.0;

}

