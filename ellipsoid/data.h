#ifndef _DATA
#define _DATA

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include "MersenneTwister.h"

#define RANDOM ( r_gen.rand() ) // RNG uniform in [0,1]
#define GRANDOM ( r_gen.randNorm(0., 1.) ) // RNG gaussian with mean 0. and variance 1.

#define PI 3.14159265358

// simulation data
#define DIMENSION 3
#define CONNECTIVITY 9.0
#define FIELD_NB 100
#define RATIO 1.5;
//#define SAMPLE_NB  10000 // number of discretization points for BP integral, per dimension: 
// 200 is ok for d=3, z=6
// 10000 is ok for d=3, z=5
#define ITERATION_NB  100 // number of iterations for the population dynamics
#define LAMBDA 0.


using namespace std;


extern MTRand r_gen;
extern ofstream data_out;

extern int fnb;
extern double *h;
extern double **n;
extern double **rxn;
extern double **r;
extern int *labels;
extern double *cumulated_weight;
extern int z;

extern int d;
extern long int samples_lim;

#endif
