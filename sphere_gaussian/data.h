#ifndef _DATA
#define _DATA

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <math.h>
#include <iostream>
#include "MersenneTwister.h"

#define RANDOM ( r_gen.rand() ) // RNG uniform in [0,1]
#define GRANDOM ( r_gen.randNorm(0., 1.) ) // RNG gaussian with mean 0. and variance 1.

#define PI 3.14159265358

// simulation data
#define DIMENSION 3
#define CONNECTIVITY 5.0
#define NORM_CUTOFF 0.01
#define BINS 1000
#define RANGE 4.
#define FIELD_NB 9
#define NA 1
#define ITERATION_NB  100 // number of iterations for the population dynamics
#define LAMBDA 0.

#define _zeroT  // comment if finite temperature

using namespace std;


class Field;

extern MTRand r_gen;
extern ofstream data_out;

extern Field *psi;
extern int fnb;
extern int bins_nb;

#endif
