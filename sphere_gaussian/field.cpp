#include <iostream>
#include "field.h"
#include "data.h"


using namespace std;


Field::Field(){
  bin_nb=BINS;
  range=RANGE;
  x_min=0.;
  x_max=x_min+range;
  dx=range/bin_nb;
 
  field=new double [bin_nb];
  double width=0.3;  
  double noise_amp=0.;
  for(int i=0;i<bin_nb;i++){
    double x=(i*range)/bin_nb;
    field[i]=(1.+noise_amp*RANDOM)*x*exp(-x*x/width/width);
  }
  set_weight();
  normalize();


  d=DIMENSION;
  n=new double [d];
  for(int i=0;i<d;i++){
    n[i]=0.;
  }
  overlap_limit=cos(PI/3.);
  updated_since_last_print=true;
  set_input_correlation(0.);
  input_fields_labels=NULL;
  input_fields_nb=0;

}
          
Field::~Field(){

  delete [] n;
  delete [] field;
  if(input_fields_labels!=NULL)
    delete [] input_fields_labels;

}


double 
Field::value(double x){
  int bin=(int)((x-x_min)*bin_nb/range);
  if(bin<0){
    cout << " Negative force? " << x << endl; exit(1);
  }
  if(bin>=bin_nb){
    return 0.;
  }
  else{
    return field[bin];
  }
}   


double
Field::norm(){
  double nnorm=0.;
  for(int i=0;i<bin_nb;i++){
    nnorm+=field[i]*dx;
  }
  return nnorm;
}
  
void
Field::normalize(){
  double nnorm=norm();
  int imax=0;
  for(int i=0;i<bin_nb;i++){
    field[i]/=nnorm;
    if(field[i]>1e-3){
      imax=i;
    }
  }
  x_max_non_zero=(double)(imax+1)*range/bin_nb;
  updated_since_last_print=true;
}

void
Field::set_weight(){
  double nnorm=norm();
  _weight=nnorm;
}

void
Field::normalize_verbose(){  // debugging
  double nnorm=0.;

   for(int i=0;i<bin_nb;i++){
    nnorm+=field[i]*dx;
    cout << field[i] << endl;
  }
   cout <<"norm " << nnorm << endl;
  double x=x_min;
  for(int i=0;i<bin_nb;i++){
    field[i]/=nnorm;
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

    out << " incoming vectors correlation " << input_correlation()<< endl;
    out << " unit vectors used for BP iteration:"  << endl;

    for(int u=0;u<d;u++)
      out << "0. ";

    for(int u=0;u<d;u++)
      out << vec(u) << " ";
    out << endl;

    for(int i=0;i<input_fields_nb;i++){
      for(int u=0;u<d;u++)
	out << "0. ";
      
      for(int u=0;u<d;u++)
	out << (psi[input_fields_labels[i]]).vec(u) << " ";
      out << endl;
    }
    
    out << endl;
    out << " Fields value: " << endl;
    
    double x=x_min+dx/10.; //avoids threshold effect due to binning
    for(int s=0; s<bin_nb;s++){
      out << x << " " << value(x) << " ";
      for(int i=0;i<input_fields_nb;i++)
	out << (psi[input_fields_labels[i]]).value(x) << " ";
      out << endl;
      x+=dx;
    }
    cout << fname<< " , weight " << weight() << endl;
    out.close();
    updated_since_last_print=false;
  }
}

void 
Field::set_unit_vector(){  // gives x_1,...,x_d coordinates of a random unit vector 
  
  if(d==3){
    double u, phi;
    
    u=2.*(RANDOM-0.5);
    phi=2.*PI*(RANDOM);
  
    n[0]=sqrt(1-u*u)*cos(phi);
    n[1]=sqrt(1-u*u)*sin(phi);
    n[2]=u;
  }
  
  if(d==2){
    double theta;
    theta=2.*PI*(RANDOM);
    n[0]=cos(theta);
    n[1]=sin(theta);
  }
}

void 
Field::print_unit_vector(){
  
  for(int u=0;u<d;u++){
    cout << 0. <<" " ;
  }

  for(int u=0;u<d;u++){
    cout << n[u] <<" " ;
  }
  cout << endl;
}

double
Field::scalar_prod(double *n2){  // tells if a sphere with unit vector n2 is in overlap with sphere with n

  double scalar=0.;

  for(int u=0;u<d;u++)
    scalar+=n[u]*n2[u];
  
  return scalar;
}


bool
Field::overlap(double *n2){  // tells if a sphere with unit vector n2 is in overlap with sphere with n

  double scalar=scalar_prod(n2);
  
  if(scalar<overlap_limit)
    return false;
  else
    return true;

}

void
Field::set_input_fields_labels(int *lab, int z){

  if(input_fields_labels!=NULL){
    delete [] input_fields_labels;
  }
  input_fields_labels= new int [z];
  input_fields_nb=z;
  
  for(int i=0;i<z;i++){
    input_fields_labels[i]=lab[i+1];
  }
}

void
Field::load_from_file(string name){
  ifstream in_stream(name.c_str());
  string s;
  int count=0;
  cout << "Read from file " << name << endl;
  while(!in_stream.eof()&& count< 100){
    getline(in_stream, s);
    if(!s.compare(" Fields value: ")){
      break;
    }
  }

  double x=x_min+dx/10.;
  double val;
  while(x<x_max){
    in_stream.ignore(256, ' ');
    in_stream >> val;
    in_stream.ignore(1024, '\n');
    set(x, val);
    x+=dx;
  }
  
  in_stream.close();
  normalize();
}
