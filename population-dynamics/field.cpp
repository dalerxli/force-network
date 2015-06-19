#include <iostream>
#include "field.h"
#include "data.h"


using namespace std;


Field::Field(){
  bin_nb=BINS;
  range=RANGE;
  x_min=MIN_FORCE;
  x_max=x_min+range;
  dx=range/bin_nb;

  field=new double [bin_nb];
  double cutoff=0.3;  
  for(int i=0;i<(int)(cutoff*bin_nb);i++){
    field[i]=RANDOM;
  }
  for(int i=(int)(cutoff*bin_nb);i<bin_nb;i++){
    field[i]=0.;
  }
  normalize();
  updated_since_last_print=true;
}

Field::~Field(){
  delete [] field;
}


double 
Field::value(double x){
  int bin=(int)((x-x_min)*bin_nb/range);
  if(bin<0){
    cout << " Negative force? " << x << endl; exit(1);
  }
  if(bin>=bin_nb){
    //    cout << " Force value ( " << x << " ) out of range ( " <<range << " )"<<  endl;
    return 0.;
  }
  else{
    return field[bin];
  }
}   
  
bool
Field::normalize(){
  double norm=0.;
  bool weird=false;
  for(int i=0;i<bin_nb;i++){
    norm+=field[i]*dx;
    if(i>0&&i<bin_nb-1){
      if(field[i]==0&&field[i+1]!=0&&field[i-1]!=0)
	weird=true;
    }
  }

  if(weird)
    cout << " weird one norm " <<  norm<< endl;
  else
    cout << " usual " << norm << endl;

  //  if(norm==0.||norm<0.1*last_norm){
  if(norm==0.||norm<0.1){
    return false;
  }
  for(int i=0;i<bin_nb;i++){
    field[i]/=norm;
  }
  updated_since_last_print=true;
  last_norm=norm;
  return true;
}

void
Field::normalize_verbose(){
  double norm=0.;

   for(int i=0;i<bin_nb;i++){
    norm+=field[i]*dx;
    cout << field[i] << endl;
  }
   cout <<"norm " << norm << endl;
  double x=x_min;
  for(int i=0;i<bin_nb;i++){
    field[i]/=norm;
    cout << field[i] << endl;
    x+=dx;  
  }
  updated_since_last_print=true;
}

void
Field::print(char* fname){
  if(updated_since_last_print){
    ofstream out;
    out.open(fname);
    double x=x_min+dx/10.; //avoids threshold effect due to binning
    for(int i=0; i<bin_nb;i++){
      out << x << " " << value(x) << endl;
      x+=dx;
    }
    cout << fname<< endl;
    out.close();
    updated_since_last_print=false;
  }
}
