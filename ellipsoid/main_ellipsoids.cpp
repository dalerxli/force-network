#include "data.h"
#include "disorder_ellipsoids.cpp"
using namespace std;

extern "C" void   dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info );
extern "C" void   dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info );
extern "C" double   dlange_(char *norm, int *m, int *n, double *a, int *lda, double *work);
extern "C" void  dgecon_(char *norm, int *n, double *a, int *lda, double *anorm, double *rcond, double *work, int *iwork, int *info);

MTRand r_gen;
ofstream data_out;
int fnb=FIELD_NB;
double *h;
double **n;
double **rxn;  // cross product r x n
double **r;
int *labels;
double *cumulated_weight;
int z;

long int samples_lim;

int d=DIMENSION;
int msize=2*DIMENSION-1;
double  U [(2*DIMENSION-1)*(2*DIMENSION-1)];
double  invU [(2*DIMENSION-1)*(2*DIMENSION-1)];

void  init_no_input(){
  h=new double [fnb];
  
  for(int i=0; i<fnb;i++){
    h[i]=RANDOM;
  }
  n=NULL;
}

void fill_U(){
  for(int u=0; u<d;u++){
    for(int v=0; v<d;v++){
      int a=(2*d-1)*u+v;
      U[a]=0.;
      for(int i=1; i<=z;i++){
	U[a]+=n[i][u]*n[i][v]/h[i];
      }
    }
    for(int v=d; v<2*d-1;v++){
      int a=(2*d-1)*u+v;
      U[a]=0.;
      for(int i=1; i<=z;i++){
	U[a]+=n[i][u]*rxn[i][v-d]/h[i];
      }
    }
  }

  for(int u=d; u<2*d-1;u++){
    for(int v=0; v<d;v++){
      int a=(2*d-1)*u+v;
      U[a]=0.;
      for(int i=1; i<=z;i++){
	U[a]+=rxn[i][u-d]*n[i][v]/h[i];
      }
    }
    for(int v=d; v<2*d-1;v++){
      int a=(2*d-1)*u+v;
      U[a]=0.;
      for(int i=1; i<=z;i++){
	U[a]+=rxn[i][u-d]*rxn[i][v-d]/h[i];
      }
    }
  }

}

void invert_U(){
  int info;
  int ipiv [2*DIMENSION-1];
  int lwork=2*DIMENSION-1;
  double work [2*DIMENSION-1];
  double work2 [4*(2*DIMENSION-1)];
  int iwork [2*DIMENSION-1];



  char *norm=new char [1];
  norm[0]='1';
  double rcond;

  double n1=dlange_(norm, &msize, &msize, U,  &msize, work); // 1-norm of U

  for(int i=0; i<(2*d-1)*(2*d-1); i++){
    invU[i]=U[i];
  }
  dgetrf_(&msize, &msize, invU, &msize, ipiv, &info); // LU factorization
  dgecon_(norm, &msize, invU, &msize, &n1, &rcond, work2, iwork, &info); // reciprocal of the condition number of U

  if(rcond<1e-10){
    cout << " non invertible matrix " << endl;
    print_vectors(); getchar();
  }

  dgetri_(&msize, invU, &msize, ipiv, work, &lwork, &info); // invert U, based on LU decomposition



  // cout << "U" << endl;
  // for(int i=0; i<2*d-1; i++){
  //     for(int j=0; j<2*d-1; j++){
  // 	cout << U[(2*d-1)*i+j] << " ";
  //     }
  //     cout << endl;
  // }
  // cout << endl;

  //  cout << "U*U^-1 " << endl;
  // for(int i=0; i<2*d-1; i++){
  //     for(int j=0; j<2*d-1; j++){
  // 	double a=0.;
  // 	for(int k=0; k<2*d-1; k++){
  // 	  a+=U[(2*d-1)*i+k]*invU[(2*d-1)*k+j];
  // 	}
  // 	cout << a << " " ;
  //     }
  //     cout << endl;
  // }

  // cout << endl << endl;

  delete [] norm;

}


void iteration(int nfl){
  free_vectors();
  set_connectivity();


  // chose the incoming/outgoing fields and there unit vectors
  set_labels(nfl);

  free_vectors();set_all_vectors();
  //  print_vectors();


  fill_U();
  invert_U();
  

  h[nfl]=0.;
  for(int u=0;u<d;u++){
    for(int v=0;v<d;v++){
      h[nfl]+=n[0][u]*invU[(2*d-1)*u+v]*n[0][v];
    }
    for(int v=d;v<2*d-1;v++){
      h[nfl]+=n[0][u]*invU[(2*d-1)*u+v]*rxn[0][v-d];
    }
  }
  for(int u=d;u<2*d-1;u++){
    for(int v=0;v<d;v++){
      h[nfl]+=rxn[0][u-d]*invU[(2*d-1)*u+v]*n[0][v];
    }
    for(int v=d;v<2*d-1;v++){
      h[nfl]+=rxn[0][u-d]*invU[(2*d-1)*u+v]*rxn[0][v-d];
    }
  }

}

void print_population(){
  data_out.seekp(ios_base::beg);
  for(int i=0;i<FIELD_NB;i++){
    data_out << i << " " << h[i] << endl;
  }
}

int main(int argc, char *argv[]){
  if ( argc != 1 ){
    fprintf(stderr, "Use with: %s \n", argv[0]);
    exit(1);
  }
  if(argc==1){
    init_no_input();
  }

  cumulated_weight= new double [fnb];
  char DATA_OUT [256];
  sprintf(DATA_OUT, "data_z%f", CONNECTIVITY);


  for(int it=0; it<ITERATION_NB*FIELD_NB; it++){
    if(it%(10*FIELD_NB)==0){
      cout << " Iterations : " << it/FIELD_NB << endl;
    }


    if(it%(10*FIELD_NB)==0){
      data_out.open(DATA_OUT);
      print_population();
      data_out.close();
    }
    int field_label=(int)(RANDOM*fnb);
    iteration(field_label);
  }


}
