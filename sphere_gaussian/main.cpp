#include <sstream>
#include "data.h"
#include "field.h"
#include "cavity.h"

using namespace std;

MTRand r_gen;
Field *psi;
BP_recursion *bp;
ofstream data_out;
int fnb;
int bins_nb;

void init_no_input(){
  fnb=FIELD_NB;
  bins_nb=BINS;
  psi=new Field [fnb];
  
  for(int i=0;i<fnb;i++){
    psi[i].normalize();
    char fname [256];
    sprintf(fname, "field_d%i_z%f_%i", DIMENSION, CONNECTIVITY, i);
    psi[i].print(fname);
  }
}

void init_with_input(string kernel, int field_nb){
  fnb=field_nb;
  bins_nb=BINS;

  psi=new Field [fnb];
  string in_file_name;
  ostringstream label;
  for(int i=0;i<fnb;i++){
    label.seekp(ios_base::beg);
    label << i;
    in_file_name=kernel+"_"+label.str();
    psi[i].load_from_file(in_file_name);
  }
  for(int i=0;i<fnb;i++){
    psi[i].normalize();
    char fname [256];
    sprintf(fname, "field_d%i_z%f_%i", DIMENSION, CONNECTIVITY, i);
    psi[i].print(fname);
  }
}

int main(int argc, char *argv[]){
  if ( argc != 1 && argc != 3 ){
    fprintf(stderr, "Use with: \n %s\n or %s\t fields_name_kernel fields_number\n", argv[0], argv[0]);
    exit(1);
  }
  if(argc==1){
    init_no_input();
  }
  else{
    init_with_input(argv[1], atoi(argv[2]));
  }

  char DATA_OUT [256];
  sprintf(DATA_OUT, "data_d_%i_z%f", DIMENSION, CONNECTIVITY);
  data_out.open(DATA_OUT);

  bp=new BP_recursion;

  for(int it=0; it<ITERATION_NB*FIELD_NB; it++){
    if(it%FIELD_NB==0){
      cout << " Iterations : " << it/FIELD_NB << endl;
    }
    int field_label=(int)(RANDOM*fnb);
    bp->iterate(field_label);

  }


  delete [] psi;
}
