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
	double cinit[nn];
	double vel[nn]; //wind reactor
	double mscl[nn];
	double flxarr[nn];
	// we don't need, right?
	double saflux[nn];
	double fc1[nn];
	double fc2[nn];

	for(i=0; i<nn; i++){
		con[i] = 1.0; //Concentration (micrograms / m3)
		cinit[i] = con[i];
		vel[i] = 0.5; //fine by Ben -> later it should vary accross the dimension, can vary magnitude and sign
		mscl[i] = 1; //map scale vector, just leave it as 1
		//everthying below this is output
		flxarr[i] = 0; //that is output from
		//saflux[i] = 0; //most meaningless but output
		//fc1[i] = 0; //even more meaningless
		//fc2[i] = 0;
	}
	for (i=0; i<5000; i++){
		//printf("[%d] Call hadvppm ... \n",i);
		//note that hadvppm( &con[0] ) is the same as hadvppm( con )
		hadvppm(nn, dx, dt, con, vel, mscl, flxarr, &flux1, &flux2);
	}

	//flxarr[nn] = 0.0;

	// save the data
	FILE * file01 = fopen("stepi.txt", "w");
	FILE * file02 = fopen("cinit.txt", "w");
	FILE * file03 = fopen("conend.txt", "w");
	FILE * file04 = fopen("flxarr.txt", "w");
	for (i=0; i<nn; i++){
		fprintf(file01, "%d\n", i);
		fprintf(file02, "%.4f\n", cinit[i]);
		fprintf(file03, "%.4f\n", con[i]);
		fprintf(file04, "%.4f\n", flxarr[i]);
	}
	fclose(file01);
	fclose(file02);
	fclose(file03);
	fclose(file04);

	return 0;
}

//this is a fucked up keyboard
/*
        Input arguments:
        nn                  Number of cells
        dt                  Time step (s)
        dx                  Length of cell (m)
        step 				dx/dt
        con                 Concentration*area vector (umol/m)
        vel                 Wind speed vector (m/s)
        mscl                Map-scale factor (squared) at cell centroid vector
        flxarr              Interfacial mass flux (umol/s) vector
        flux1-2             Boundary fluxes (umol/s)
                             (1=west/south, 2=east/north)
        saflux              Interfacial mass flux times area vector
                             (used for tracer transport)
        fc1-2               Conc*area vector change from flux (umol/m)
                             (1=west/south, 2=east/north)
                             (used for Process Analysis)

		Output arguments:
        con
        flxarr
        flux1-2
		saflux

*/


void hadvppm(HADVPPM_ARGS){
    //
    int i;
    double TWO3RDS = 2.0/3.0;
    double Conp;
    double Conm;
    double x;

    // BUFFER MATRICES
    double fm[nn];
    double fp[nn];
    //double fc1[nn];
    //double fc2[nn];
    //double saflux[nn];
    double cm[nn];
    double cl[nn];
    double cr[nn];
    double dc[nn];
    double c6[nn];

    for (i=0; i<nn; i++){
        fm[i] = 0;
        fp[i] = 0;
        cm[i] = 0;
        cl[i] = 0;
        cr[i] = 0;
        dc[i] = 0;
        c6[i] = 0;
        //fc1[i] = 0;
        //fc2[i] = 0;
    }
    // Be caurfull: The matrices in C start by 0 while in Fortran by 1
    cm[1]  = con[1];
    cm[nn-1] = con[nn-2];
    cm[2] = (con[2] + con[1])/2;
    cm[nn-2] = (con[nn-2] + con[nn-3])/2;

    for(i=2; i<nn-2; i++){
        dc[i] = 0.5*(con[i+1] - con[i-1]);
        Conp = con[i+1] - con[i];
        Conm = con[i] - con[i-1];
        if ((Conp)*(Conm) > 0){
            dc[i] = copysign(1,dc[i]) * fmin(abs(dc[i]),fmin(2.0*abs(Conp),2.0*abs(Conm)));
        } else {
			dc[i] = 0;
        }
        /*if (i != 3)
			cm[i] = con[i-1] + 0.5*(con[i] - con[i-1]) + (dc[i-1] - dc[i])/6.0;
  */}

    for(i=2; i<nn-3; i++){ // ToDo: Merge this for-loop with the previous
		cm[i+1] = con[i] + 0.5*(con[i+1] - con[i]) + (dc[i] - dc[i+1])/6.0;
    }

    for(i=1; i<nn-1; i++){
        cr[i] = cm[i+1];
        cl[i] = cm[i];
    }

    for(i=1; i<nn-1; i++){
    	if ((cr[i] - con[i])*(con[i] - cl[i])>0) {
			dc[i] = cr[i] - cl[i];
			c6[i] = 6.0*(con[i] - 0.5*(cl[i] + cr[i]));
			if (dc[i]*c6[i] > dc[i]*dc[i]){
				cl[i] = 3.0*con[i] - 2.*cr[i];
			} else if (-dc[i]*dc[i] > dc[i]*c6[i]){
				cr[i] = 3.0*con[i] - 2.*cl[i];
			}
    	} else {
			cl[i] = con[i];
			cr[i] = con[i];
    	}
        dc[i] = cr[i] - cl[i];
        c6[i] = 6.0*(con[i] - 0.5*(cl[i] + cr[i]));
    }

    for(i=1; i<nn-1; i++){
		x = fmax(0,-vel[i-1]*(dt/dx));
		//printf("[%d]I AM x1: %f\n",i,x);
		fm[i] = x*(cl[i] + 0.5*x*(dc[i] + c6[i]*(1.0 - TWO3RDS*x)));
        //if (x >= 1)
			//printf("Courant number %f is bigger than 1", x);
		x = fmax(0,vel[i]*(dt/dx)); // ToDo: maybe this is not necessary
		//printf("[%d]I AM x2: %f\n",i,x);
		//if (x >= 1)
			//printf("Courant number %f is bigger than 1", x);
        fp[i] = x*(cr[i] - 0.5*x*(dc[i] - c6[i]*(1.0 - TWO3RDS*x)));
	}

	if (vel[0] > 0){
		x = vel[0]*(dt/dx);
        fp[0] = x*con[0];
	}

	if (vel[nn-2] < 0){
		x = -vel[nn-2]*(dt/dx);
        fp[nn-1] = x*con[nn-1];
	}

	flxarr[0] = (fp[0] - fm[1])*(dx/dt);
	//printf("flxarr:%f, fp0:%f, fm1:%f \n", flxarr[0], fp[0], fm[1]);

	//saflux[0] = flxarr[0]*step;

	for (i=1; i<nn-1; i++){
		//printf("[%d]", i);
		flxarr[i] = (fp[i] - fm[i+1])*(dx/dt);
		con[i] = con[i] - mscl[i]*(flxarr[i] - flxarr[i-1])*(dt/dx);
	}

	*flux1 = mscl[1]*flxarr[0];
    *flux2 = mscl[nn-2]*flxarr[nn-2];

	/* test: print output arguments
    printf("mscl\n");
    for (i = 0; i<nn; i++){
		printf("%.4f\n", mscl[i]);
    }
    printf("the Flux1:%f\n", flux1);
    printf("the Flux2:%f\n", flux2);
    printf("flxarr		con		   fm		fp\n");
    for (i = 0; i<nn; i++){
		printf("%.4f		%.4f		%.4f		%.4f\n", flxarr[i],  con[i], fm[i], fp[i]);
    }
    printf("cm		cl		cr		dc		c6\n");
    for (i = 0; i<nn; i++){
		printf("%.4f		%.4f		%.4f		%.4f		%.4f\n", cm[i], cl[i], cr[i], dc[i], c6[i]);
    }
    printf("\n END\n");
    */
}

/* We don't need it .....

void mnratio(MNRATIO_ARGS){
	int ii;
	int i;
	int spec;
	int j;
	int k;
	int flag; // 1 => x-advection, 2 => y-advection
	int order[nspc];
	double c1d[nn];
	double c2d[nn];
	double mn[nn];
	double eps = 1.0e-20;

	if (ii>=ngas){
		spec = (ii-ngas)%naero;
		if (spec == 1){
			int l = 0;
			for (i=i1; i<i2; i++){
				l = l + 1;
				c2d[l] = c1d[l];
				if (c2d[l] <= 0){
					if (c2d[l] > -eps){
						c2d[l] = eps;
					} else {
						printf("In mnratio, number concentrations are less than -eps\n");
						printf("Species = %f\n", spname(order(ii));
						printf("In mnratio, number concentrations are less than -eps");
					}
				}
			}
		}
	}

}
*/


