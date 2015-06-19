#include "data.h"

using namespace std;

void free_vectors(){

  if(n!=NULL){
    for(int i=0;i<z+1;i++){
      delete [] n[i];
    }
    delete [] n;
  }
  n=NULL;

  if(rxn!=NULL){
    for(int i=0;i<z+1;i++){
      delete [] rxn[i];
    }
    delete [] rxn;
  }
  rxn=NULL;

  if(r!=NULL){
    for(int i=0;i<z+1;i++){
      delete [] r[i];
    }
    delete [] r;
  }
  r=NULL;
}

void ellipsoid_RNG(double *rand_vec, double minor_axis_length, double major_axis_length){ // random vectors on a d-dimensional ellipsoid, with d-1 minor axis, and 1 major axis

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
  
bool is_label_ok(int lb){
  
  for(int i=0;i<lb;i++){
    if(labels[i]==labels[lb])
      return false;
  }
  return true;
}

void set_unit_vector(double *nn, double *rr, double *rxnn){  // gives x_1,...,x_d coordinates of a random unit vector 

  double minor_axis=1.;  // minor=1.93; major=1. is M&M's from Chaikin
  double major_axis=RATIO;
  ellipsoid_RNG(rr, minor_axis, major_axis);


  double norm=0.;
  for(int i=0; i<d-1;i++){
    nn[i]=rr[i]/minor_axis/minor_axis;
    norm+=nn[i]*nn[i];
  }    
  nn[d-1]=rr[d-1]/major_axis/major_axis;
  norm+=nn[d-1]*nn[d-1];

  norm=sqrt(norm);

  for(int i=0; i<d;i++){
    nn[i]/=norm;
  }

  if(d==2){
    rxnn[0]=rr[0]*nn[1]-rr[1]*nn[0];
  }
  if(d==3){
    rxnn[0]=rr[1]*nn[2]-rr[2]*nn[1];
    rxnn[1]=-rr[0]*nn[2]+rr[2]*nn[0];
  }


  if(d!=2&&d!=3){
    cout << "set_unit_vector(double *nn, double *rxnn) not implemented for d="<< d << endl; 
    exit(1);
  }
    
}

double scalar_prod(double *n1, double *n2){

  double scalar=0.;

  for(int u=0;u<d;u++)
    scalar+=n1[u]*n2[u];
  
  return scalar;
}


bool overlap(double *n1, double *n2){  // tells if a sphere with unit vector n2 is in overlap with sphere with n ---> to be adapted for ellipsoids.

  double overlap_limit=cos(PI/3.);
  double scalar=scalar_prod(n1, n2);
  
  if(scalar<overlap_limit)
    return false;
  else
    return true;

}

bool is_vector_ok(int lb){
  
  for(int i=0;i<lb;i++){
    if(overlap(n[lb], n[i]))
      return false;
  }

  return true;
}


void set_labels(int nfl){

  // computing cumulated weights. Here is stupid as they all have the same weight, but can be useful if we want to give different weights
  cumulated_weight[0]=1.; 
  for(int i=1;i<fnb;i++){
    cumulated_weight[i]=cumulated_weight[i-1]+1.;
  }

  for(int i=0;i<fnb;i++){
    cumulated_weight[i]/=cumulated_weight[fnb-1];
  }


  if(labels!=NULL)
    delete [] labels;
  labels=new int [z+1];  // the list of incoming fields
  labels[0]=nfl;
  for(int i=1;i<=z;i++){
    do{
      double r=RANDOM;
      for(int j=0;j<fnb;j++){
	if(cumulated_weight[j]>r){
	  labels[i]=j;
	  break;
	}
      }

    }while(!is_label_ok(i));
  }

}



bool geometrical_constraint_ok(){  // this tests if there is a hemisphere free of incoming vectors.

  double orth_vector [d];

  if(d==2){
    for(int i=0;i<=z;i++){

      orth_vector[0]=-n[i][1];
      orth_vector[1]=n[i][0];
      
      bool has_pos=false;
      bool has_neg=false;
      double cutoff=1e-7;  // should be zero. if positive, we are asking that there is no free cap with angle arccos(cutoff).

      for(int j=0;j<=z;j++){
	if(scalar_prod(orth_vector, n[j])<-cutoff)
	  has_neg=true;
	if(scalar_prod(orth_vector, n[j])>cutoff)
	  has_pos=true;
      }
      if(!has_pos||!has_neg)
	return false;
    }    
    return true;
  }


  if(d==3){
    for(int i=0;i<=z;i++){
      for(int j=i+1;j<=z;j++){

	orth_vector[0]=n[i][1]*n[j][2]-n[i][2]*n[j][1];
	orth_vector[1]=-n[i][0]*n[j][2]+n[i][2]*n[j][0];
      	orth_vector[2]=n[i][0]*n[j][1]-n[i][1]*n[j][0];

	double norm=0.;	
	for(int u=0;u<d;u++){
	  norm+=orth_vector[u]*orth_vector[u];
	}
	norm=sqrt(norm);
	for(int u=0;u<d;u++){
	  orth_vector[u]/=norm;
	}


	bool has_pos=false;
	bool has_neg=false;
	double cutoff=1e-7;  // should be zero. if positive, we are asking that there is no free cap with angle arccos(cutoff).
	for(int k=0;k<=z;k++){
	  if(scalar_prod(orth_vector, n[k])<-cutoff)
	    has_neg=true;
	  if(scalar_prod(orth_vector, n[k])>cutoff)
	    has_pos=true;
	}
	if(!has_pos||!has_neg)
	  return false;
      }
    }
    return true;
  }
  cout << "test of vector set is not implemented for d="<< d << endl;
  exit(-1);
  return false;
}

void set_outgoing_vector(){

  do{
    set_unit_vector(n[0], r[0], rxn[0]);
  }while(!is_vector_ok(z)||!geometrical_constraint_ok());

}



void set_all_vectors(){ // sets incoming and outgoing vectors

  n=new double* [z+1];
  for(int i=0;i<=z;i++){
    n[i]=new double [d];
  }
  r=new double* [z+1];
  for(int i=0;i<=z;i++){
   r[i]=new double [d];
  }
  rxn=new double* [z+1];
  for(int i=0;i<=z;i++){
    rxn[i]=new double [d-1];
  }


  do{
    for(int i=0;i<=z;i++){
      do{
	set_unit_vector(n[i], r[i], rxn[i]);
      }while(!is_vector_ok(i));
    }

  }while(!(geometrical_constraint_ok()));


}

void print_vectors(){

  cout << endl << "  unit vectors :" << endl;
  for(int i=0;i<=z;i++){
    for(int u=0; u<d; u++){
      cout << 0. << " ";
    }
      
    for(int u=0; u<d; u++){
      cout << n[i][u] << " ";
    }
    for(int u=0; u<d; u++){
      cout << r[i][u] << " ";
    }

    for(int u=0; u<d-1; u++){
      cout << rxn[i][u] << " ";
    }


    cout << endl;
  }
}


void set_connectivity(){
    z=CONNECTIVITY-1;
  /*
  if(d==2){// stupid thing: P(z)=0.5 for z_iso, and adjust P(z-1) and P(z+1) to have the correct connectivity: valid only for coneectivities b/n 3.5 and 4.5
    double p_of_z[3];
    p_of_z[1]=0.5;
    p_of_z[2]=0.5*(CONNECTIVITY-3.5);
    p_of_z[0]=0.5-p_of_z[2];

    double s=RANDOM;
    double cumul=0.;

    for(int i=0;i<3; i++){
      cumul+=p_of_z[i];
      if(s<cumul){
	z=i+2;
	if(z==3){
	  samples_lim=10000;
	}
	if(z==4){
	  samples_lim=40000;
	}
	if(z==5){
	  samples_lim=200000;
	}
	break;
      }
    }
  }

  if(d==3){//  gaussian with cutoffs at min_z and min_z+range
    int range=5;
    int min_z=7;

    double * p_of_z=new double [range]; // from 4 to 11
    double norm=0.;
    double var=4.;

    for(int i=0; i<range;i++){
      p_of_z[i]=exp(-(i+min_z-CONNECTIVITY)*(i+min_z-CONNECTIVITY)/var);
      norm+=p_of_z[i];
    }
    p_of_z[0]/=norm;
    for(int i=1; i<range;i++){
      p_of_z[i]=p_of_z[i-1]+p_of_z[i]/norm;
    }
    
    double choice=RANDOM;
    for(int i=0; i<range;i++){
      if(choice<p_of_z[i]){
	z=i+min_z;
	break;
      }
    }
  }
  */
}
