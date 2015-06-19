#ifndef _DATA
#define _DATA

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <math.h>
#include <iostream>
#include "MersenneTwister.h"

#define RANDOM ( r_gen.rand() )
#define PI 3.14159265358

// simulation data
#define DIMENSION 2 
#define BINS 200
#define RANGE 1.
#define MIN_FORCE 0.
#define FIELD_NB 40
#define SAMPLE_NB  1000 // number of discretization points for BP integral
#define ITERATION_NB  10000 // number of iterations for the population dynamics
#define BETA 1.
#define LAMBDA 4.  // =beta*lambda if _zeroT

#define _zeroT  // comment if finite temperature

using namespace std;


class Field;

extern MTRand r_gen;
extern double beta;
extern double lambda;


extern Field *psi;
#endif
