/*****
 *
 *  header.h
 *
 **/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.c"

//arguments to function hadvppm
#define HADVPPM_ARGS int nn, double step, double con[], double vel[], double mscl[], double flxarr[], double* flux1, double* flux2, double saflux[], double fc1[], double fc2[]
#define MNRATIO_ARGS // TODO
#define DEBUG 0

// DO NOT CHANGE
#define MX1D 97
#define PI 3.141592653589793238462643383

// Prototype
void hadvppm(HADVPPM_ARGS);
