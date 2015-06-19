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
  double range;
  double dx;
  bool updated_since_last_print;
  double last_norm;
 public:

  double value(double);

  void set(double x, double a){
    int bin=(int)((x-x_min)*bin_nb/range);
    field[bin]=a;
    updated_since_last_print=true;
  }

  Field();
  ~Field();

  bool normalize();
  void normalize_verbose();
  void print(char*);
  
};


#endif

