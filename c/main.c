/* this is the equivalent to xyadvec.f */

#include "header.h"

int main(){

	int nn = 90;
	int i;
	double dt = 900;
	double dx = 36000;
	double flux1 = 0; //fluxes at the boundaries output from the subroutine
	double flux2 = 0;
	double con[nn];
	double vel[nn]; //wind reactor
	double mscl[nn];
	double flxarr[nn];
	double saflux[nn];
	double fc1[nn];
	double fc2[nn];

	for(i=0; i<nn; i++){
		con[i] = sin(2.0/(double)nn * PI * (double)i);
		vel[i] = 0.5; //fine by Ben -> later it should vary accross the dimension, can vary magnitude and sign
		mscl[i] = 1; //map scale vector, just leave it as 1
		//everthying below this is output
		flxarr[i] = 0; //that is output from
		saflux[i] = 0; //most meaningless but output
		fc1[i] = 0; //even more meaningless
		fc2[i] = 0;
	}
	//note that hadvppm( &con[0] ) is the same as hadvppm( con )
    hadvppm(nn, dx, dt, con, vel, mscl, flxarr, &flux1, &flux2, saflux, fc1, fc2);

	return 0;
}
