//#include<openmp.h> //is this the correct compiler directive?

//this is a fucked up keyboard

//arguments to function hadvppm
#define HADVPPM_ARGS int nn, double dt, double dx, double con[], double vel[], double mscl[], double flxarr[], double* flux1, double* flux2, double saflux[], double fc1[], double fc2[]
/*
        Input arguments:
        nn                  Number of cells
        dt                  Time step (s)
        dx                  Length of cell (m)
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
*/
//function prototype
void hadvppm(HADVPPM_ARGS){
    //
    int i;
    double TWO3RDS = 2.0/3.0;
    double Conp;
    double Conm;
    double x;
    double step = dt/dx; // ToDo

    //DUMMY VARIABLES

    int MX1D = 1;

    // BUFFER MATRICES
    double fm[MX1D];
    double fp[MX1D];
    double cm[MX1D];
    double cl[MX1D];
    double cr[MX1D];
    double dc[MX1D];
    double c6[MX1D];

    for (i=0; i<nn; i++){ // Why nn?? The size of matrices fm, fp are MX1D
        fm[i] = 0;
        fp[i] = 0;
        fc1[i] = 0;
        fc2[i] = 0;
    }
    // Be caurfull: The matrices in C start by 0 while in Fortran by 1
    cm[1]  = con[1];
    cm[nn-1] = con[nn-2];
    cm[2] = (con[2] + con[1])/2;
    cm[nn-2] = (con[nn-2] + con[nn-3])/2;

    for(i=3; i<nn-2; i++){
        dc[i] = 0.5*(con[i+1] - con[i-1]);
        Conp = con[i+1] - con[i];
        Conm = con[i] - con[i-1];
        if ((Conp)*(Conm) > 0){
            dc[i] = copysign(1,dc[i]) * fmin(abs(dc[i]),fmin(2.0*abs(Conp),2.0*abs(Conm))); // sign ?????
        } else {
			dc[i] = 0;
        }
    }

    for(i=3; i<nn-3; i++){ // ToDo: Merge this for-loop with the previous
		cm[i+1] = con[i] + 0.5*(con[i+1] - con[i]) + (dc[i] - dc[i+1])/6.0;
    }

    for(i=2; i<nn-1; i++){
        cr[i] = cm[i+1];
        cl[i] = cm[i];
    }

    for(i=2; i<nn-1; i++){
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

    for(i=2;i<nn-1;i++){
		x = fmax(0,-vel[i-1]*(dt/dx));
		fm[i] = x*(cl[i] + 0.5*x*(dc[i] + c6[i]*(1.0 - TWO3RDS*x)));
        if (x >= 1)
			printf("Courant number %f is bigger than 1", x);
		x = fmax(0,vel[i-1]*(dt/dx)); // ToDo: maybe this is not necessary
		if (x >= 1)
			printf("Courant number %f is bigger than 1", x);
        fp[i] = x*(cr[i] - 0.5*x*(dc[i] - c6[i]*(1.0 - TWO3RDS*x)));
	}

	if (vel[1] > 0){
		x = vel[1]*(dt/dx);
        fp[1] = x*con[1];
	}

	if (vel[nn-1] > 0){
		x = -vel[nn-1]*(dt/dx);
        fp[nn] = x*con[nn];
	}

	flxarr[1] = (fp[1] - fm[2])*(dx/dt);
	saflux[1] = flxarr[1]*(dt/dx);

	for (i=1; i<nn-1; i++){
		flxarr[i] = (fp[i] - fm[i+1])*(dx/dt);
		con[i] = con[i] - mscl[i]*(flxarr[i] - flxarr[i-1])*(dt/dx);
        saflux[i] = flxarr[i]*(dt/dx);
        fc1[i] =   mscl[i]*flxarr[i-1]*(dt/dx);
        fc2[i] = - mscl[i]*flxarr[i]*(dt/dx);
	}

	*flux1 = mscl[2]*flxarr[1];
    *flux2 = mscl[nn-1]*flxarr[nn-1];

    printf("hello");
}

