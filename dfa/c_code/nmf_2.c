#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <time.h>

#define M_PI 3.14159265358979323846
#define PI 3.14159265
#define LEN(arr) ((sizeof (arr) / sizeof (*arr))) //macro to define the size of an array
#define tol 0.0001


float drand()   /* uniform distribution, (0..1] */
{
  return (rand()+1.0)/(RAND_MAX+1.0);
}

float random_normal()
 /* normal distribution, centered on 0, std dev 1 */
{
  return sqrt(-2*log(drand())) * cos(2*M_PI*drand());
}

//generate a random matrix
float** rand_matrix( float **mat, int m, int n, float mult_elem){
  //m is running throug the rows and n through the column
  //mult_elem is the element we want to moltiply with
    float** backup = mat;
    for (int i=0;i<m;i++){
        for (int j=0;j<n;j++){
            mat[i][j] = mult_elem* (1.0 + 0.5*random_normal());
        }
    }
    return backup;
}

void print_matrix( float **mat, int row, int col){
    printf("\nMatrix dimensions: [%d,%d]\n", row, col);

    for (int i=0;i<row;i++){
        for (int j=0;j<col;j++){
            printf(" %f ", mat[i][j]);
        }
        printf("\n");
    }
}

//Ht  X_columns x n_components
//H   n_components x X_columns
///row=X_Columns col=n_components
float** transpose_matrix( float **mat, float **matT, int row, int col){
    float** backup = matT;
    float** backup2 = mat;

    for (int i =0; i < col; i++){
      for (int j=0; j< row; j++){
        backup[j][i]=backup2[i][j];
      }
    }

    return backup;

}
//typedef float fftw_complex[2];
float define_avg_element(float **X, int X_rows, int X_columns, int n_components){
  //X_rows --> n_samples
  //X_columns --> n_compoments
  float summa = 0.0 ;
  int n_elems = X_rows*X_columns;
  for (int i=0; i< X_rows; i++){
    for (int j=0; j<X_columns;j++){
      //printf("%.2f\n", X[i][j]);
      summa+=X[i][j];
    }
  }
  if (summa<0.0){
    summa = -1*summa;
  }
  printf("Final sum %.4f\n",summa);
  float avg = sqrt( (summa/n_elems)/ n_components) ;

  printf("%.2f\n",avg);
  return avg;
}


float 

float* nmf(int n_components,
            int max_iter,
            int random_state,
            float *X,
            int X_rows,
            int X_columns,
            float *W_cython,
            float *H_cython,
            int verbose)
{

  //transform X into a matrix
  float **Xmat;
  if (verbose==1){
    printf("Conversion of Cython vector signal X to matrix");
  }
  Xmat = (float **) malloc(X_rows * sizeof(float*));
  for (int i =0; i< X_rows; i++){
    Xmat[i] = (float *)malloc(X_columns*sizeof(float));
  }
  int counter =0 ;
  for (int i=0; i<X_rows;i++){
    for(int j=0; j<X_columns;j++){
      Xmat[i][j]=X[counter];
      counter+=1;
    }
  }


  //take the average elmeent
  float avg_elem;
  avg_elem = define_avg_element(Xmat, X_rows, X_columns, n_components);
  //initialize random H
  float **H;
  H = (float **) malloc(n_components * sizeof(float*)) ;
  for(int i = 0; i<n_components; i++) {
        H[i] = (float *) malloc(X_columns * sizeof(float));
  }
  if (verbose==1){
    printf("Initializing matrix H random\n");
  }
  H  = initialize_H(H, avg_elem, X_rows, X_columns, n_components, random_state);
  printf("Starting H Matrix");
  print_matrix(H, n_components, X_columns );
  if (verbose==1){
    printf("Initialized H matrix\n");
    //print_matrix(H, n_components, X_columns );
  }
  //initialize random W
  float **W;
  W = (float **)malloc(X_rows *sizeof(float*));
  for (int i=0; i<X_rows; i++){
      //X_rows = n_samples   X_columns = n_features
      W[i]=(float *) malloc(n_components *sizeof(float));
  }
  W = initialize_W(W, avg_elem, X_rows, X_columns, n_components, random_state);
  printf("Starting W matrix");
  print_matrix(W, X_rows, n_components);
  if (verbose==1){
    printf("Initialized W matrix\n");
    //print_matrix(W, X_rows, n_components);
    printf("Coordinate descent\n");
  }

  //TODO regularization, skip for the moment
  //Coordinate descent solver
  coordinate_descent_solver(H, W, Xmat, n_components, X_rows, X_columns,max_iter, verbose);
  if (verbose==1){
    printf("W matrix\n");
    print_matrix(W,X_rows, n_components);
    printf("H matrix\n");
    print_matrix(H, n_components, X_columns);
  }
  //printf("Final W matrix\n");
  //print_matrix(W,X_rows, n_components);
  //printf("Final H matrix\n");
  //print_matrix(H, n_components, X_columns);
  //vectorize the result to be sent to cython
  int W_counter = 0;
  for (int i =0; i< X_rows; i++){
    for (int j=0; j<n_components; j++){
      W_cython[W_counter] = W[i][j];
      W_counter+=1;
    }
  }

  int H_counter =0;
  for (int i =0; i< n_components; i++){
    for (int j=0; j<X_columns; j++){
      //printf("%.2f\n",H[i][j]);
      H_cython[H_counter] = H[i][j];
      H_counter+=1;
    }
  }
  //return W_cython, H_cython;
  free(H);
  free(W);
  return 0;
}
