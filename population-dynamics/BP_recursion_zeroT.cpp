#include "BP_recursion.h"
#include "data.h"
#include "field.h"

#ifdef _zeroT

void 
BP_recursion::random_vector(double* rvec, int d){  // gives x,y,z coordinates of a random unit vector 
  
  if(d==3){
    double u, phi;
    
    u=2.*(RANDOM-0.5);
    phi=2.*PI*(RANDOM);
  
    rvec[0]=sqrt(1-u*u)*cos(phi);
    rvec[1]=sqrt(1-u*u)*sin(phi);
    rvec[2]=u;
  }
  
  if(d==2){
    double theta;
    theta=2.*PI*(RANDOM);
    rvec[0]=cos(theta);
    rvec[1]=sin(theta);
  }
}


bool 
BP_recursion::generate_forces(double x, double y [], double* new_n, double **n, int d, int z){

  for(int i=0;i<z;i++){
    y[i]=MIN_FORCE+RANDOM*RANGE;
  }
  
  // impose mechanical equilibrium
  // requires to select d particles among neighbors, and solve the system of d equations required by mechanical equilibrium. Then we can set the forces associated with these d particles

  int ran_neighbor [d];
  bool is_selected [z];
  for(int i=0;i<z;i++){
    is_selected[i]=false;
  }

  for(int i=0;i<d;i++){
    do{
	do{
	    ran_neighbor[i]=(int)(RANDOM*z);
	}while(ran_neighbor[i]==z);
		
	is_selected[ran_neighbor[i]]=false;
	for(int j=0;j<i;j++){
	    if(ran_neighbor[i]==ran_neighbor[j]) is_selected[ran_neighbor[i]]=true;
	}
    }while(is_selected[ran_neighbor[i]]);
    is_selected[ran_neighbor[i]]=true;
  }
  
  // now we find the mechanical equilibrium

  double sum_of_forces [d];
  for(int u=0;u<d;u++){
    sum_of_forces[u]=x*new_n[u];
    for(int i=0;i<z;i++){
      if(!is_selected[i])
	sum_of_forces[u]+=y[i]*n[i][u];
    }
  }
  
  if(d==2){
    int a=ran_neighbor[0];
    int b=ran_neighbor[1];
    y[b]=(-sum_of_forces[1]+sum_of_forces[0]*n[a][1]/n[a][0]);
    y[b]/=(n[b][1]-n[b][0]*n[a][1]/n[a][0]);
    
    y[a]=(-sum_of_forces[0]-y[b]*n[b][0])/n[a][0];
  }
  else{
    cout << "Mechanical equilibrium for d=" << d << " not implemented yet" << endl; exit(1);
  }
  for(int i=0;i<d;i++){
//    if(y[ran_neighbor[i]]<0.||y[ran_neighbor[i]]>(MIN_FORCE+RANGE))
    if(y[ran_neighbor[i]]<0.)
      return false;
  }

  return true;
}

void
BP_recursion::print_state(double x, double y [], double* new_n, double **n, int z,  int d){
  ofstream out("state.dat");
  out << endl << endl <<" state : "<< endl;
  
  for(int u=0;u<d;u++){
      out << 0. << " ";
  }
  for(int u=0;u<d;u++){
    out << x*new_n[u] << " ";
  }
  out <<endl;
  for(int i=0;i<z;i++){
    for(int u=0;u<d;u++){
      out << 0. << " ";
    }
    for(int u=0;u<d;u++){
      out << y[i]*n[i][u] << " ";
    }
    out << endl;
  }

  double sum_of_forces [d];
  for(int u=0;u<d;u++){
    sum_of_forces[u]=x*new_n[u];
    for(int i=0;i<z;i++){
      sum_of_forces[u]+=y[i]*n[i][u];
    }
  }
  
  for(int u=0;u<d;u++){
    out << 0. << " ";
  }
  for(int u=0;u<d;u++){
    out << sum_of_forces[u] << " ";
  }
  out <<endl;

  out.close();
}

double 
BP_recursion::BP_integral(double x, double* new_n, double **n, int z, int labels [], int d){

  double y [z];   // forces in the selected fields
  double integral=0.;
  int samples=0;

  while(samples<SAMPLE_NB){
    double part_integral=0.;
    int count=0;

    while(!generate_forces(x, y, new_n, n, d, z)){
      count++;
      if(count>1000){
	  print_state(x, y, new_n, n, z, d); return 0.; // we could not find any mechanical equilibrium between chosen particles
      }
      //      if(count%100==0)
      //	cout << " sample " << samples << " trying to find equilibrium " << count << endl;
    }

    double sum1=x*x;
    for(int i=0;i<z;i++){
      sum1+=y[i]*y[i];
    }
    part_integral=exp(-lambda*sum1);

    for(int i=0;i<z;i++){
      part_integral*=psi[labels[i]].value(y[i]);
    }
    integral+=part_integral;
    samples++;
  }
  return integral/samples;

}
  
void 
BP_recursion::assign_vectors(double* new_n, double** n, int z, int d){

  for(int i=0;i<z;i++){
    n[i]=new double [d];
    random_vector(n[i], d);  // assign a random vector to the selected incoming fields
  }
  
  bool accept=false;
  do{
    random_vector(new_n, d);
    double scal=0.;
    for(int i=0;i<z;i++){ // ensure that at least one of them is in the hemisphere opposite to new_n
      for(int u=0;u<d;u++){ 
	scal=n[i][u]*new_n[u];
      }
      if(scal<0.)
	accept=true;
    }
  }while(!accept);

}

void 
BP_recursion::iterate(int new_field_label){


  int d=DIMENSION;

  int count=0;
  do{
    int z=3; // connectivity-1
    int labels [z];  // the list of incoming fields
    for(int i=0;i<z;i++){
      do{
	labels[i]=(int)(RANDOM*FIELD_NB);
      }while(labels[i]==new_field_label);
    }
    
    
    double** n;  // unit vectors of the incoming fields
    n=new double* [z];
    double* new_n;  // unit vector of the new field
    new_n=new double [z];
    
    int angular_samples=0;
    do{
	assign_vectors(new_n, n, z, d); // average over angles too!    
	
	double x=MIN_FORCE;
	double dx=RANGE/BINS;
	x+=dx/10.; //avoids threshold effect due to binning
	while(x<MIN_FORCE+RANGE){
	    if(angular_samples==0)
		(psi[new_field_label]).set(x, BP_integral(x, new_n, n, z, labels, d));
	    else
		(psi[new_field_label]).set(x, (psi[new_field_label]).value(x)+BP_integral(x, new_n, n, z, labels, d));
	    x+=dx;
	}
	angular_samples++;
    }while(angular_samples<100);

    count++;
  }while(!(psi[new_field_label]).normalize());  // !psi.normalize() means that the new field is identically zero! Happens sometimes if neighbors are not well chosen

  char fname [256];
  sprintf(fname, "fields/field_beta%lf_lamda%lf_%i",beta, lambda, new_field_label);
  (psi[new_field_label]).print(fname);
}

#endif
