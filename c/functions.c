#include<stdio.h>
#include<math.h>
//#include<openmp.h> //is this the correct compiler directive?

//this is a fucked up keyboard

//arguments to function hadvppm
#define HADVPPM_ARGS int nn, double dt, double dx, double* con[], double vel[], double mscl[], double* flxarr[], double* flux1, double* flux2, double* saflux[], double* fc1[], double* fc2[]
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
    double TWO3RGS = 2/3;
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
        if ((con[i+1] - con[i])*(con[i] - con[i-1]) > 0){
                dc[i] =
        }
    }

    printf("hello");
}

int main(){
    return 0;
}

