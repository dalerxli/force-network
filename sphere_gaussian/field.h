#ifndef _FIELD
#define _FIELD
#include <fstream>
#include "data.h"

using namespace std;

class Field
{
 private:

  
  double* field;
  int bin_nb;
  double x_min;
  double x_max;
  double x_max_non_zero;
  double range;
  double dx;
  bool updated_since_last_print;
  double _weight;
  
  int d; //dimension
  double *n; // unit vector associated with the field

  double overlap_limit;
  
  double input_correlation_value;

  int * input_fields_labels;
  int input_fields_nb;
  

 public:

  
  double value(double);

  void set(double x, double a){
    int bin=(int)((x-x_min)*bin_nb/range);
    field[bin]=a;
    updated_since_last_print=true;
  }


  void reset(){
    for(int bin=0; bin<bin_nb; bin++)
      field[bin]=0.;

    updated_since_last_print=true;
  }
  
  double max_non_zero(){
    return x_max_non_zero;
  }
  
  double weight(){
    return _weight;
  }
  void set_weight();

  Field();
  ~Field();

  void load_from_file(string);
  double vec(int u){
    return n[u];
  }

  void set_unit_vector();
  void print_unit_vector();
  
  bool overlap(double*);
  double scalar_prod(double*);
  double norm();
  void normalize();
  void normalize_verbose();
  void print(char*);

  void set_input_correlation(double icv){
    input_correlation_value=icv;
  }
  double input_correlation(){
    return input_correlation_value;
  }

  void set_input_fields_labels(int *, int);
};


#endif

