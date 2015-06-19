#include "data.h"
#include "field.h"
#include "BP_recursion.h"

using namespace std;

MTRand r_gen;
double beta=BETA;
double lambda=LAMBDA;
Field *psi;

BP_recursion *bp;

int main(int argc, char *argv[]){
  if (argc !=1 ){
    fprintf(stderr, "Use with: %s\n", argv[0]);
    exit(1);
  }

  bp=new BP_recursion;
  psi=new Field [FIELD_NB];
  for(int i=0;i<FIELD_NB;i++){
    psi[i].normalize();
    char fname [256];
    sprintf(fname, "fields/field_beta%lf_lamda%lf_%i",beta, lambda, i);
    psi[i].print(fname);
  }

  for(int it=0; it<ITERATION_NB; it++){
    if(it%10==0){
      cout << " Iterations : " << it << endl;
    }
    int field_label=(int)(RANDOM*FIELD_NB);
    bp->iterate(field_label);

  }


  delete [] psi;
}
