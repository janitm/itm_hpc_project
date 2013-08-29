/* this is the equivalent to xyadvec.f */

#include "header.h"

int main(){

	int nn = 30;
	int i;
	double dt = 0.1;
	double dx = 0.1;
	double flux1 = 1;
	double flux2 = 1;
	double con[nn];
	double vel[nn];
	double mscl[nn];
	double flxarr[nn];
	double saflux[nn];
	double fc1[nn];
	double fc2[nn];

	for(i=0; i<nn; i++){
		con[i] = sin(2.0/(double)nn * PI * (double)i);
		vel[i] = 0.5;
		mscl[i] = 1;
		flxarr[i] = 1;
		saflux[i] = 1;
		fc1[i] = 1;
		fc2[i] = 1;

		//note that hadvppm( &con[0] ) is the same as hadvppm( con )
        hadvppm(nn, dt, dx, con, vel, mscl, flxarr, &flux1, &flux2, saflux, fc1, fc2);
        }

	return 0;
}

//#define HADVPPM_ARGS int nn, double dt, double dx, double con[], double vel[], double mscl[], double flxarr[], double* flux1, double* flux2, double saflux[], double fc1[], double fc2[]
