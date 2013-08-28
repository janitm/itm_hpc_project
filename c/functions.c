#include<stdio.h>
#include<math.h>
//#include<openmp.h> //is this the correct compiler directive?


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
    double TWO3RGS = 2/3;
    
    printf("hello");
}

int main(){
    return 0;
}

