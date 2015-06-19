#ifndef _BPR
#define _BPR

#include <fstream>

#include<gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

#include "data.h"

using namespace std;

class BP_recursion{


 private:
  
  ofstream out_test;
  void BP_integral(double*, double*);
  void generate_forces(double**, gsl_vector**, double*, double**, bool*, bool*);
  void print_state(double []);
  bool is_label_ok(int);
  bool is_vector_ok(int);
  bool geometrical_constraint_ok();
  void init_labels(int);
  void set_all_vectors();
  void set_outgoing_vector();
  void lin_solve(double [], int [], bool []);

  double interaction_entropy;
  double site_entropy;
  double edge_entropy;
  long int iterations;
  long int samples_lim;

  void analyse(double, double); // debugging function
  int z;
  int *labels;
  double **n;
  void record_vectors();
  void print_vectors();
  void free_vectors();
  int d;
  double set_connectivity();

  bool fail_flag;
  void compute_entropy(int);

  void unitball_RNG(double*, int);
  double * y0;
  double rad_0;
  gsl_matrix * q_matrix;
  gsl_matrix * invtr_matrix;
  void compute_algebra();
  void allocate_algebra();
  void free_algebra();

  
  double * cumulated_weight;
 public:
  void iterate(int);
  BP_recursion(){
    labels=NULL;
    n=NULL;
    interaction_entropy=0.;
    site_entropy=0.;
    edge_entropy=0.;
    iterations=0;
    compute_entropy(0);
    out_test.open("forces.dat");
    cumulated_weight=new double [fnb];
  }

  ~BP_recursion(){
    delete [] labels;
    delete [] n;
    delete [] cumulated_weight;
  }

};





#endif

