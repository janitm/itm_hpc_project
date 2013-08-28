#include<stdio.h>
#include<math.h>
//#include<openmp.h> //is this the correct compiler directive?


//arguments to function hadvppm
#define HADVPPM_ARGS double* nn, float dt, float dx, double* con, double* vel, double mscl,double* flxarr, double flux1,double flux2,double saflux,double fc1,double fc2

//function prototype
void hadvppm(HADVPPM_ARGS)
{
    printf("hello");
}

int main()
{
    return 0;
}

