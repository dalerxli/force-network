#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "MersenneTwister.h"

#define RANDOM ( r_gen.rand() ) // RNG uniform in [0,1]
#define GRANDOM ( r_gen.randNorm(0., 1.) ) // RNG gaussian with mean 0. and variance 1.

using namespace std;


MTRand r_gen;
int d=2;


void ellipsoid_RNG(double *rand_vec, double minor_axis_length, double major_axis_length){

  double q=0.;
  double r=0.;
  double u;
  do{
    q=0.;
    r=0.;
  
    
    for(int i=0;i<d-1;i++){
      rand_vec[i]=GRANDOM;
      q+=rand_vec[i]*rand_vec[i];
      r+=rand_vec[i]*rand_vec[i]/minor_axis_length/minor_axis_length;
    }
    rand_vec[d-1]=GRANDOM;
    q+=rand_vec[d-1]*rand_vec[d-1];
    r+=rand_vec[d-1]*rand_vec[d-1]/major_axis_length/major_axis_length;
  
    q=sqrt(q);
    r=sqrt(r);
    
    u=RANDOM;
    
  }while(minor_axis_length*r < q*u);

  for(int i=0;i<d-1;i++){
    rand_vec[i]*=minor_axis_length/q;
  }
  rand_vec[d-1]*=major_axis_length/q;

}


int main(int argc, char *argv[]){
  if ( argc != 1 ){
    fprintf(stderr, "Use with: %s \n", argv[0]);
    exit(1);
  }
  
  double *nn;
  nn=new double [d];

  while(1){
    ellipsoid_RNG(nn, 1., 3.);
    
    for(int u=0; u<d; u++){
      cout << nn[u] << " ";
    }
    cout << endl;
  }
  delete [] nn;
}
