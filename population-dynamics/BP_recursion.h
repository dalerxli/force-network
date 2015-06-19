#ifndef _BPR
#define _BPR

using namespace std;

class BP_recursion{


 private:
  
  void random_vector(double*, int);
  double BP_integral(double, double*, double**, int , int [], int);
  void assign_vectors(double*, double**, int, int);
  bool generate_forces(double, double [], double*, double**, int, int);
  void print_state(double, double [], double*, double**, int , int);
 public:
  void iterate(int);
  BP_recursion(){}
  ~BP_recursion(){}


};
#endif
