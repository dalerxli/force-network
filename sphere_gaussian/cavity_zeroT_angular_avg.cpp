#include <stdlib.h>
#include "cavity.h"
#include "data.h"
#include "field.h"


#ifdef _zeroT


void
BP_recursion::compute_algebra(){

  // this function computes the matrices necessary to get the force balance.
  // we need to find solutions to N*f_in=-f_out, where N is a d*z matrix called incoming_matrix in the code
  // the solutions of this system are obtained by a QR decomposition of the transpose t_incoming_matrix of N.
  // As each time we change f_out, the solution changes of course, BUT the incoming matrix depends only 
  // on the incoming directions n. This function thus gives only the two matrices Q and (R_1^T)^-1, where 
  // R_1 is the d*d non-zero part of the R matrix. From these two matrices
  // it is straightforward to get f_in knowing f_out:
  // f_in = Q y, with y = [ - (R_1^T)^-1 f_out, y2 ]
  // y2 being an arbitrary vector of dimension z-d
 
  // fill the matrix N^T (called t_incoming_matrix)
  gsl_matrix * t_incoming_matrix = gsl_matrix_calloc (z, d);
  for(int i=0;i<z;i++){
    for(int u=0;u<d;u++){
      gsl_matrix_set(t_incoming_matrix, i, u, (psi[labels[i+1]]).vec(u));
    }
  }


  gsl_matrix * r_large_matrix=  gsl_matrix_alloc(z, d); // the matrix R
  gsl_matrix * r_small_matrix=  gsl_matrix_alloc(d, d);  // the non-zero part of R
  gsl_vector * tau = gsl_vector_alloc(d);   // a vector that GSL uses for the QR decomposition

  gsl_linalg_QR_decomp(t_incoming_matrix, tau); gsl_linalg_QR_unpack (t_incoming_matrix, tau, q_matrix, r_large_matrix);


  // now we have the Q and R matrices of the QR decomposition of the transpose of the incoming matrix
  // we now need to extract the non-zero part of R
  for(int i=0;i<d;i++){
    for(int u=0;u<d;u++){
      double a=gsl_matrix_get(r_large_matrix, i, u);
      gsl_matrix_set(r_small_matrix, i, u, a);
    }
  }

  // now we invert r_small_matrix, using LU decomposition: it is easy because r_small_matrix is already an upper triangular matrix. The decomposition is then trivial: r_small_matrix=I*r_small_matrix

  gsl_permutation * perm=gsl_permutation_alloc (d); // the LU decomposition contains a permutation P (A is such that PA=LU)
  gsl_permutation_init(perm); // of course, here P is identity
  gsl_linalg_LU_invert(r_small_matrix, perm, invtr_matrix);

  // and at last we transpose, to get the result
  gsl_matrix_transpose(invtr_matrix);

  // the last thing we need to compute now is the coordinates of the center of the cube [0, f_max]^z, in the y-basis
  // (ie the basis related to the original forces basis via y = Q^T * forces )
  // this point, called y0 here, is used to generate the random points in the y-basis for the integration
  gsl_vector * Y0 = gsl_vector_alloc(z);
  gsl_vector * F0 = gsl_vector_alloc(z);
  for(int i=0;i<z; i++){
    gsl_vector_set(F0, i, 0.5*(psi[labels[i+1]]).max_non_zero());
  }
  gsl_blas_dgemv(CblasTrans, 1., q_matrix, F0, 0., Y0);
  for(int i=0;i<z; i++){
    y0[i]=gsl_vector_get(Y0, i);
  }

  // desallocate
  gsl_vector_free(tau);
  gsl_vector_free(F0);
  gsl_vector_free(Y0);
  gsl_matrix_free(t_incoming_matrix);
  gsl_matrix_free(r_large_matrix);
  gsl_matrix_free(r_small_matrix);
  gsl_permutation_free(perm);

}


void
BP_recursion::print_state(double y []){
  ofstream out("state.dat");

  out << endl << endl <<" forces : "<< endl;
  for(int i=0;i<=z;i++){
    out << y[i] << endl;
  }
  
  out << endl << endl <<" state : "<< endl;
  
  for(int i=0;i<=z;i++){
     for(int u=0;u<d;u++){
       out << 0. << " ";
     }
    for(int u=0;u<d;u++){
      out << y[i]*(psi[labels[i]]).vec(u) << " ";
    }
    out << endl;
  }

  double sum_of_forces [d];
  for(int u=0;u<d;u++){
    sum_of_forces[u]=0.;
      for(int i=0;i<=z;i++){
      sum_of_forces[u]+=y[i]*(psi[labels[i]]).vec(u);
    }
  }
  out << "sum of forces ";
   for(int u=0;u<d;u++){
     out << 0. << " ";
   }
  for(int u=0;u<d;u++){
    out << sum_of_forces[u] << " ";
  }
  out <<endl;

   out.close();
}
 

void 
BP_recursion::free_vectors(){


  if(n!=NULL){
    for(int i=0;i<z+1;i++){
      delete [] n[i];
    }
    delete [] n;
  }
  n=NULL;
}

void 
BP_recursion::record_vectors(){


  n=new double* [z+1];  
  for(int i=0;i<=z;i++){
    n[i]=new double [d];
    for(int u=0;u<d;u++){
      n[i][u]=(psi[labels[i]]).vec(u);
    }
  }

}  
void
BP_recursion::unitball_RNG(double *rand_vec, int dim){
  double norm=0.;
  for(int i=0;i<dim;i++){
    rand_vec[i]=GRANDOM;
    norm+=rand_vec[i]*rand_vec[i];
  }
  norm=sqrt(norm);
  
  double u=pow(RANDOM, 1./dim);
  for(int i=0;i<dim;i++){
    rand_vec[i]*=u/norm;
  }
  
}


void
BP_recursion::generate_forces(double ** f, gsl_vector ** y_inc_full, double * rad_E_0, double ** y_E_0, bool * non_neg_forces, bool *out_of_range){

  /*********************/
  // note: the nice thing with having Q (q_matrix) orthogonal is that uniform sampling over the (z-d)-dimensional 
  // affine plane of allowed forces (that we note S) in the y-basis GUARANTEES uniform sampling of this plane in the f-basis.
  // 
  // Now, as we want the forces to belong to [0,fmax], it is enough to sample the intersection U of [0, fmax]^z with the affine plane 
  // S (which is defined by Q(y_inc_part, y_tail), y_tail in [reals]^(z-d)).
  // In the y-basis, this means sampling the image "Q^T(U)" of U through Q^T. 
  // Unfortunately, it is difficult in practice to compute Q^T(U).
  //
  // A simpler (although non optimal) strategy, is to sample a finite space E that we know contains Q^T(U).
  // A way to do this, is to take E as being the intersection of S with the ball containing [0, fmax]^z (in the f-basis), which is 
  // centered on y0 and has radius rad0.
  // E is then a ball of dimension (z-d), and we can easily compute its radius and center. We obtain the center position in the
  // y-basis, and we then generate random points in E, using a standard unit-ball uniform RNG.
  // This is the strategy we adopt here.
  //
  // It should work fine if (z-d) is small, ie when the ratio of the volume of U on the volume of E is reasonnable. For large z-d,
  // we will lose a lot of time generating random points outside the interesting intersection U
  /********************/
  


  //  we generate random vector on E, and put it in the z-d last values of y_inc_full
  double rand_vec [z-d];
  unitball_RNG(rand_vec, z-d);
  gsl_vector * forces;

  forces= gsl_vector_alloc(z);
  for(int i=0;i<BINS;i++){
    if(!out_of_range[i]){
      for(int j=d;j<z;j++){
	gsl_vector_set(y_inc_full[i], j, y_E_0[i][j]+rad_E_0[i]*rand_vec[j-d]);
      }
      // now we set y_inc_full to be the set of forces that satisfies force balance through forces = q_matrix * y_inc_full
      

      gsl_blas_dgemv (CblasNoTrans, 1., q_matrix, y_inc_full[i], 0., forces);
      
      non_neg_forces[i]=true;

      double *fi=f[i];
      for(int j=1;j<=z;j++){
	fi[j]=gsl_vector_get(forces, j-1);
	// if(i<2)
	//   cout << "forces for i="<< i << " : " << fi[j] << endl;
	if(fi[j]<0.)
	  non_neg_forces[i]=false;
      }
    }
  }
  gsl_vector_free(forces);  

}


void
BP_recursion::BP_integral(double *x, double *field_value){


  fail_flag=false;
  double **f;   // forces in the selected fields

  f=new double* [BINS];

  /**** BEGINNING OF THE ALGEBRA PART *****/
  gsl_vector ** y_inc_part;
  gsl_vector ** y_inc_full;
  y_inc_part=new gsl_vector* [BINS];
  y_inc_full=new gsl_vector* [BINS];

  double ** y_E_0;
  y_E_0=new double* [BINS];
  double *dist2_y0_y_E_0;
  dist2_y0_y_E_0=new double [BINS];
  double *rad_E_0;
  rad_E_0=new double [BINS];
  bool *out_of_range;
  out_of_range=new bool [BINS];
  int nb_in_range=0;

  for(int i=0; i<BINS;i++){
    f[i]=new double [z+1];
    y_E_0[i]=new double [z];

    for(int j=1;j<=z;j++){
      f[i][j]=0.;
    }
    f[i][0]=x[i];

    y_inc_part[i] = gsl_vector_alloc(d);
    y_inc_full[i] = gsl_vector_alloc(z);

  // set y_inc_part to be (- f_outgoing)
    for(int u=0;u<d;u++){
      gsl_vector_set(y_inc_part[i], u, -f[i][0]*n[0][u]);
    }


  // y_inc_part = invtr_matrix * (- f_outgoing), as required by QR decomposition
    gsl_blas_dtrmv (CblasLower, CblasNoTrans, CblasNonUnit, invtr_matrix, y_inc_part[i]);


  // fill y_inc_full with y_inc_part for the d first values
    for(int u=0;u<d;u++){
      gsl_vector_set(y_inc_full[i], u, gsl_vector_get(y_inc_part[i], u));
    }

  // determine the center y_E_0, and the radius rad_E_0 of the sphere E (see the note on function generate_forces)
  // note: from an implementation point of view, it is surely inefficient to declare and fill y_E_0, 
  // as its coordiantes are simple and easy to access, but I keep it at the moment for the clarity of the algebra.

    for(int j=0;j<d;j++){
      y_E_0[i][j]=gsl_vector_get(y_inc_full[i], j);
    }
    for(int j=d;j<z;j++){
      y_E_0[i][j]=y0[j];
    }

    dist2_y0_y_E_0[i]=0.;
    for(int u=0;u<d;u++){
      dist2_y0_y_E_0[i]+=(y0[u]-y_E_0[i][u])*(y0[u]-y_E_0[i][u]);
    }



    out_of_range[i]=false;
    //    if(i<2)
      //      cout << sqrt(dist2_y0_y_E_0[i]) << " " << rad_0 << endl;
    if(sqrt(dist2_y0_y_E_0[i])>rad_0){ // there is no solution to sum f=0 in [0, f_max]^z, so field_value(x[i]) will be set to 0, using out_of_range
      out_of_range[i]=true;
      rad_E_0[i]=0.;
    }
    else{
      rad_E_0[i]=sqrt(rad_0*rad_0-dist2_y0_y_E_0[i]);
      nb_in_range++;
    }

  }
  
  if(nb_in_range==1){  // means that only one binning for f_out (maybe) works: this is a tricky case, that we avoid: by setting all out_of_range[i]=true, the new field will be identically zero, and we will generate a new set of unit vector
    for(int i=0; i<BINS;i++){
      out_of_range[i]=true;
    }
    cout << " tricky case, skipping ... " << endl;
    nb_in_range=0;
  }
  

  /**** END OF THE ALGEBRA PART *****/


  double integral [BINS];
  double part_integral [BINS];
  bool * non_neg_forces;
  non_neg_forces= new bool [BINS];

  for(int i=0; i<BINS;i++){
    integral[i]=0.;
  }

  long int samples=0;

  //  cout << nb_in_range << endl;

  while(samples<samples_lim&&nb_in_range>0){
    //    cout << samples << endl;

    generate_forces(f, y_inc_full, rad_E_0, y_E_0, non_neg_forces, out_of_range);
    
    for(int i=0; i<BINS;i++){

      if(!out_of_range[i]&&non_neg_forces[i]){
	part_integral[i]=1.;

	for(int j=1;j<=z;j++){
	  part_integral[i]*=psi[labels[j]].value(f[i][j]);
	}
	integral[i]+=part_integral[i];
      }
    }
    samples++;
  }


  for(int i=0; i<BINS;i++){
    gsl_vector_free(y_inc_part[i]);
    gsl_vector_free(y_inc_full[i]);
    delete [] f[i];
    delete [] y_E_0[i];
  }

  for(int i=0; i<BINS;i++){
    field_value[i]=(integral[i]/samples_lim);
  }


  delete [] y_inc_part;
  delete [] y_inc_full;
  delete [] f;
  delete [] y_E_0;
  delete [] dist2_y0_y_E_0;
  delete [] rad_E_0;
  delete [] out_of_range;
  delete [] non_neg_forces;

  
}
  
bool
BP_recursion::is_label_ok(int lb){
  
  for(int i=0;i<lb;i++){
    if(labels[i]==labels[lb])
      return false;
  }
  return true;
}

bool
BP_recursion::is_vector_ok(int lb){
  
  double vector [d];

  for(int u=0;u<d;u++){
    vector[u]=(psi[labels[lb]]).vec(u);
  }
  
  for(int i=0;i<lb;i++){
    if((psi[labels[i]]).overlap(vector))
      return false;
  }
  return true;
}


void
BP_recursion::init_labels(int nfl){
  cumulated_weight[0]=(psi[0]).weight();
  for(int i=1;i<fnb;i++){
    cumulated_weight[i]=cumulated_weight[i-1]+(psi[i]).weight();
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

void
BP_recursion::set_outgoing_vector(){

  
  do{
    (psi[labels[0]]).set_unit_vector();
  }while(!is_vector_ok(z)||!geometrical_constraint_ok());

}


bool
BP_recursion::geometrical_constraint_ok(){
  double orth_vector [d];

  if(d==2){
    for(int i=0;i<=z;i++){

      orth_vector[0]=-(psi[labels[i]]).vec(1);
      orth_vector[1]=(psi[labels[i]]).vec(0);
      
      bool has_pos=false;
      bool has_neg=false;
      double cutoff=0.3;
      for(int j=0;j<=z;j++){
	if((psi[labels[j]]).scalar_prod(orth_vector)<-cutoff)
	  has_neg=true;
	if((psi[labels[j]]).scalar_prod(orth_vector)>cutoff)
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
      
	orth_vector[0]=(psi[labels[i]]).vec(1)*(psi[labels[j]]).vec(2)-(psi[labels[i]]).vec(2)*(psi[labels[j]]).vec(1);
	orth_vector[1]=-(psi[labels[i]]).vec(0)*(psi[labels[j]]).vec(2)+(psi[labels[i]]).vec(2)*(psi[labels[j]]).vec(0);
	orth_vector[2]=(psi[labels[i]]).vec(0)*(psi[labels[j]]).vec(1)-(psi[labels[i]]).vec(1)*(psi[labels[j]]).vec(0);
      
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
	double cutoff=0.;  // should be zero. if positive, we are asking that there is no free cap with angle arccos(cutoff).
	for(int k=0;k<=z;k++){
	  if((psi[labels[k]]).scalar_prod(orth_vector)<-cutoff)
	    has_neg=true;
	  if((psi[labels[k]]).scalar_prod(orth_vector)>cutoff)
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

void
BP_recursion::set_all_vectors(){ // sets incoming and outgoing vectors


  int count=0;
  bool fail=false;
  do{
    count=0;
    for(int i=0;i<=z;i++){
      do{
	(psi[labels[i]]).set_unit_vector();
	count++;
      }while(!is_vector_ok(i)&&count<1000);
    }
    if(count>=1000)
      fail=true;
      
  }while(!(geometrical_constraint_ok())&&!fail);

  double input_correl=0.;  // a correlation in incoming vectors directions

  for(int u=0; u<d; u++){
    double udir=0.;
    for(int i=1;i<=z;i++){
      udir+=(psi[labels[i]]).vec(u);
    }
    input_correl+=udir*udir;
  }

  input_correl=sqrt(input_correl);
  (psi[labels[0]]).set_input_correlation(input_correl);
  (psi[labels[0]]).set_input_fields_labels(labels, z);

}

void
BP_recursion::print_vectors(){

  cout << endl << "  unit vectors :" << endl;
  for(int i=0;i<=z;i++){
    for(int u=0; u<d; u++){
      cout << 0. << " ";
    }
      
    for(int u=0; u<d; u++){
      cout << (psi[labels[i]]).vec(u) << " ";
    }
    cout << endl;
  }
}

void
BP_recursion::allocate_algebra(){
  q_matrix = gsl_matrix_alloc(z, z);
  invtr_matrix  = gsl_matrix_alloc(d, d);
  y0=new double [z];

}
void
BP_recursion::free_algebra(){
  gsl_matrix_free(q_matrix);
  gsl_matrix_free(invtr_matrix);
  delete [] y0;
}


double
BP_recursion::set_connectivity(){
  
  if(d==2){// stupid thing: P(z)=0.5 for z_iso, and adjust P(z-1) and P(z+1) to have the correct connectivity: valid only for coneectivities b/n 3.5 and 4.5
    double p_of_z[3];
    int min_z=2;
    p_of_z[1]=0.5;
    p_of_z[2]=0.5*(CONNECTIVITY-3.5);
    p_of_z[0]=0.5-p_of_z[2];

    if(z==2){
      samples_lim=10000;
    }
    if(z==3){
      samples_lim=50000;
    }
    if(z==4){
      samples_lim=200000;
    }
    return 0.5*(z+1)*(p_of_z[z-min_z]);
  }

  if(d==3){
    int min_z=3;
    int range=5;
    double * p_of_z=new double [range]; // from 4 to 11
    double norm=0.;
    double var=1.5;

    for(int i=0; i<range;i++){
      p_of_z[i]=exp(-(i+min_z+1-CONNECTIVITY)*(i+min_z+1-CONNECTIVITY)/var);
      norm+=p_of_z[i];
    }
    /*
    p_of_z[0]/=norm;
    for(int i=1; i<range;i++){
      p_of_z[i]=p_of_z[i-1]+p_of_z[i]/norm;
    }
    */
    for(int i=0; i<range;i++){
      p_of_z[i]/=norm;
    }
    if(z==3){
      samples_lim=10000;
    }
    if(z==4){
      samples_lim=50000;
    }
    if(z==5){
      samples_lim=200000;
    }
    if(z>=6){
      samples_lim=400000;
    }
    return 0.5*(z+1)*(p_of_z[z-min_z]);
  }

  return 1.;
}

void 
BP_recursion::iterate(int nfl){


  d=DIMENSION;

  int count=0;
  Field new_psi;
  bool fail;

  double * x;
  x=new double [BINS];
  double * integral;
  integral=new double [BINS];
  double dx=RANGE/bins_nb;
  x[0]=0.1*dx; // we shift with 0.9*dx to avoid threshold effect due to binning
  for(int i=1;i<BINS;i++){
    x[i]=x[i-1]+dx;
  }
  


  z=d;  // d+1-1
  int z_max;
  if(d==2)
    z_max=4;
  if(d==3)
    z_max=7;

  int ang_samples = 0;
  int ang_samples_lim = NA;
  (psi[nfl]).reset();
  do{// z-loop
    double weight=set_connectivity();
    init_labels(nfl);
    do{// angular average loop
      do{
	// chose the incoming/outgoing fields and there unit vectors

	free_vectors();set_all_vectors();record_vectors();
	
	// set up the algebra of the QR decomposition
	allocate_algebra();
	compute_algebra();
	rad_0=0.;
	for(int i=1;i<=z;i++){
	  rad_0+=0.25*(psi[labels[i]]).max_non_zero()*(psi[labels[i]]).max_non_zero();
	}
	rad_0=sqrt(rad_0);
	
	cout << "connectivity : " << z+1 << "; angle set nb:" << ang_samples<< endl;

	
	// iterate
	BP_integral(x, integral);
	for(int i=0;i<BINS;i++){
	  new_psi.set(x[i], integral[i]);
	}
	
	
	// now test if the iteration is acceptable: 
	// - the new field has to be non-zero, 
	//  // - the new field should not show a tendancy to shrink a delta peak (this is an ad hoc test to prevent a field \delta(x) to occur if the unit vectors are not compatible with force balance
	fail=true;
	double norm=new_psi.norm();
	if(norm>NORM_CUTOFF){
	  
	  int new_max_non_zero=0;
	  for(int i=0;i<BINS;i++){
	    (psi[labels[0]]).set(x[i], (psi[labels[0]]).value(x[i])+weight*new_psi.value(x[i]));
	    // if((psi[labels[0]]).value(x[i])/norm>1e-3)
	    //   new_max_non_zero=i;
	  }
	  
	  
	  fail=false;
	  // for(int j=1; j<=z;j++){
	  //   double cutoff=0.2; // this cutoff is difficult to set: should try different values
	  //   if(x[new_max_non_zero]>cutoff*(psi[labels[j]]).max_non_zero()){	  
	  //     fail=false;
	  //   }
	  // }
	  
	  // // if(fail) at this point means that this is probably difficult to get force balance: f_out (and all f_in) can only be 0 (or almost): this does not happen in a real packing
	  // if(fail){
	  //   cout << "fail due to shrinking field, ";
	  //   (psi[labels[0]]).print("temp"); //getchar();
	  // }
	}
	else{
	  cout << "fail due to low norm ( " << norm << " ), ";
	  if(isnan(norm)){
	    new_psi.print("temp"); getchar();
	  }
	}
	
	if(fail)
	  cout << "trying ... " << endl;      
	
	
	free_algebra();
	
	
	count++;
	
      }while((psi[labels[0]]).norm()<NORM_CUTOFF||fail); // 1 angle
      ang_samples++;
    }while(ang_samples<ang_samples_lim);
    free_vectors();
    z++;
    ang_samples=0;
  }while(z<z_max);

  (psi[labels[0]]).set_weight();
  compute_entropy(labels[0]);  // don't normalize before computing the entropy

  (psi[labels[0]]).normalize();

  char fname [256];
  sprintf(fname, "field_d%i_z%f_%i",DIMENSION, CONNECTIVITY, labels[0]);
  (psi[labels[0]]).print(fname);

  delete [] x;
  delete [] integral;

}



void
BP_recursion::compute_entropy(int a_field){  // has to called BEFORE normalizing psi[a_field]

  // **** edge entropy 
  double psi_norm=(psi[a_field]).norm();
  
  double dx=RANGE/bins_nb;
  double x=dx/10.;
  while(x<RANGE){
    (psi[a_field]).set(x, (psi[a_field]).value(x)*exp(LAMBDA*x*x));  // now psi is \hat{psi}, modulo normalization
    x+=dx;
  }
  double hat_psi_norm=(psi[a_field]).norm();  

  double e_entropy=0.;
  x=dx/10.;
  while(x<RANGE){// \int dx \hat{psi}[a_field]*\hat{psi}[another_field] e^{-lanbda*x*x}
    e_entropy+=(psi[a_field]).value(x)*(psi[a_field]).value(x)*exp(-LAMBDA*x*x)*dx/hat_psi_norm/psi_norm;  
    x+=dx;
  }
  e_entropy=log(e_entropy);
  
  
  // ***** interaction entropy  
  

  double int_entropy=log(hat_psi_norm);
  
  
  // ***** site entropy  
  int another_field;
  do{
    another_field=(int)(RANDOM*fnb);
  }while(another_field==a_field);

  x=dx/10.;
  while(x<RANGE){
    (psi[another_field]).set(x, (psi[another_field]).value(x)*exp(LAMBDA*x*x)); // now psi is \hat{psi}, modulo normalization
    x+=dx;
  }
  double hat_psi_norm2=(psi[another_field]).norm();
  
  
  double s_entropy=0.;
  x=dx/10.;
  while(x<RANGE){// \int dx \hat{psi}[a_field]*\hat{psi}[another_field] e^{-lanbda*x*x}
    s_entropy+=(psi[another_field]).value(x)*(psi[a_field]).value(x)*exp(-LAMBDA*x*x)*dx/hat_psi_norm/hat_psi_norm2;  
    x+=dx;
  }
  s_entropy=log(s_entropy);
  
  iterations++;
  interaction_entropy+=int_entropy;
  site_entropy+=s_entropy;
  edge_entropy+=e_entropy;


// now backwards: \hat{psi} ----> \psi
  x=dx/10.;
  while(x<RANGE){
    (psi[another_field]).set(x, (psi[another_field]).value(x)*exp(-LAMBDA*x*x));
    (psi[a_field]).set(x, (psi[a_field]).value(x)*exp(-LAMBDA*x*x));
    x+=dx;
  }
  (psi[another_field]).normalize();
  (psi[a_field]).normalize();

  cout << "S_int: "<<  interaction_entropy/iterations << "   S_site: " << site_entropy/iterations << " S_edge " << edge_entropy/iterations << endl;
  data_out  << "S_int: "<<  interaction_entropy/iterations << "   S_site: " << site_entropy/iterations << " S_edge " << edge_entropy/iterations << endl;

}

void
BP_recursion::analyse(double x, double dx){  // debugging function
  if((psi[labels[0]]).value(x)<0.01*(psi[labels[0]]).value(x-dx)&&(psi[labels[0]]).value(x-dx)>0.1){
    cout << "stop at " << x << endl;
    cout << "z is " << z << endl;
    cout << "previous value " <<  (psi[labels[0]]).value(x-dx) << " " << (psi[labels[0]]).value(x-2.*dx) << endl;
    cout << " current value " << (psi[labels[0]]).value(x) << endl;
    cout << "incoming fields ";
    for(int i=1;i<=z;i++)
      cout << labels[i] << " ";
    cout << endl;
    for(int i=0;i<=z;i++)
      (psi[labels[i]]).print_unit_vector();

    cout << endl;
    getchar();
  }
}


#endif
