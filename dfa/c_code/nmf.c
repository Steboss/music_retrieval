#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <time.h>

#define M_PI 3.14159265358979323846
#define PI 3.14159265
#define LEN(arr) ((sizeof (arr) / sizeof (*arr))) //macro to define the size of an array
#define tol 0.000001


double drand()   /* uniform distribution, (0..1] */
{
  return (rand()+1.0)/(RAND_MAX+1.0);
}

double random_normal()
 /* normal distribution, centered on 0, std dev 1 */
{
  return sqrt(-2*log(drand())) * cos(2*M_PI*drand());
}

//generate a random matrix
double** rand_matrix( double **mat, int m, int n, double mult_elem){
  //m is running throug the rows and n through the column
  //mult_elem is the element we want to moltiply with
    double** backup = mat;
    for (int i=0;i<m;i++){
        for (int j=0;j<n;j++){
            mat[i][j] = mult_elem* (1.0 + 0.5*random_normal());
        }
    }
    return backup;
}

void print_matrix( double **mat, int row, int col){
    printf("\nMatrix dimensions: [%d,%d]\n", row, col);

    for (int i=0;i<row;i++){
        for (int j=0;j<col;j++){
            printf(" %f ", mat[i][j]);
        }
        printf("\n");
    }
}

double** transpose_matrix( double **mat, double **matT, int row, int col){
    double** backup = matT;
    double** backup2 = mat;

    for (int i =0; i < col; i++){
      for (int j=0; j< row; j++){
        backup[j][i]=backup2[i][j];
      }
    }

    return backup;

}

double ** multiply_matrix(double **mat, double **mat2, double ** result, int row, int col){

  double ** backup = result;
  for (int i =0; i< row; i++){
    for(int j =0; j< col; j++){
      for (int k=0; k< col; k++){
        backup[i][j]+= mat[i][k]*mat2[k][j];
      }
    }
  }
  return backup;
}
//typedef double fftw_complex[2];
double define_avg_element(double **X, int X_rows, int X_columns, int n_components){
  //X_rows --> n_samples
  //X_columns --> n_compoments
  double summa = 0 ;
  int n_elems = X_rows*X_columns;
  for (int i=0; i< X_rows; i++){
    for (int j=0; j<X_columns;j++){
      summa+=X[i][j];
    }
  }
  double avg = sqrt( (summa/n_elems)/ n_components) ;
  printf("%.2f\n",avg);
  return avg;
}

double** initialize_H(double** H, double avg_elem, int X_rows, int X_columns, int n_components, int r_seed){
  //default n_components = 2
  //we're assuming that X_coluimns is the number of n_features

  //generate a random matrix starting from avg
  srand(4*r_seed);
  //the matrix is n_components * n_features --> n_features is the number of columns in X
  H=rand_matrix(H, n_components, X_columns, avg_elem);
  //absolute value
  for (int i =0; i<n_components; i++){
    for (int j=0; j<X_columns;j++){
      H[i][j] = fabs(H[i][j]);
    }
  }
  //printf("H matrix\n");
  //print_matrix(H, n_components, X_columns);
  return H;
}

double** initialize_W(double ** W, double avg_elem, int X_rows, int X_columns, int n_components, int r_seed){

  //initialise W
  srand( (10*r_seed));
  //W is n_samples * n_compoments, where n_samples is the umber of rows in X
  W= rand_matrix(W, X_rows, n_components, avg_elem);
  for (int i =0; i<X_rows; i++){
    for (int j=0; j<n_components;j++){
      W[i][j] = fabs(W[i][j]);
    }
  }
  //print_matrix(W, X_rows, n_components);
  return W;

}
/* Update coordinate
X W Ht Python -->  X, HT, W in C
X X_columns * X_rows  (n_features/components * n_samples)
W n_samples * n_components
Ht n_samles * n_components

A n_components is shape[1] from HT
1. Transpose HT --> HTT
2. dot HTT and HT --> n_components * n_components
3. dot X and Ht --> X-rows *n_components
4. Permutations so


In the second case we want ot update H so  we pass
XT, HT, W --> XT W ? in C

A n_components is shape[1] from W  (so the same as above)
1. Transpose of W --> WT
2. dot WT and W
3. dot XT and W
4. Permutations

*/
double update_coordinate_descent(double** Xmat, double** Ht, double** W, int n_components, int X_rows, int X_columns){
  /*Update W to minimize the objective function, iterating once over all
    coordinates. By symmetry, to update H, one can call
    _update_coordinate_descent(X.T, Ht, W, ...)*/

  //take the transpose of Ht
  double **HtT;
  HtT = (double ** )malloc(n_components*sizeof(double*));
  for (int i=0; i<n_components; i++){
    HtT[i] = (double*)malloc(X_columns*sizeof(double));
  }
  HtT = transpose_matrix(Ht, HtT, n_components, X_columns);
  printf("HtT :\n");
  print_matrix(HtT, n_components, X_columns);

  //Create a 2D array to  store the values
  //Now the multiplication between  H and Ht will have X_columns*X_columns
  //to be safe let's put is = n_components * n_componets
  double **HttHt;
  HttHt = (double**) malloc(n_components*sizeof(double*));
  for (int i=0;i<n_components;i++){
    HttHt[i]=(double*)malloc(n_components*sizeof(double));
  }
  //initialize also the result X*Ht
  double **XHt;
  XHt = (double**) malloc(X_rows*sizeof(double*));
  for (int i=0;i<X_rows;i++){
    XHt[i]=(double*)malloc(n_components*sizeof(double));
  }


  double HttHtcurr_summ = 0.0;

  for (int i=0; i<n_components; i++){
    for(int j=0; j<n_components;j++){
      for(int k=0; k< n_components;k++){
        //TODO: Check why we have a numerical discrepancy with Python here
        HttHtcurr_summ+=HtT[i][k]*Ht[k][j];
      }
      //fill the array
      HttHt[i][j] = HttHtcurr_summ;
      HttHtcurr_summ=0;

    }

  }

  double XHtcurr_summ= 0.0; //does this hold for everything?
  //Ht will have n_components * n_compoments
  //X is X_rows * n_compoments (X_features)
  //row * column
  for (int i =0; i<X_rows; i++){
    for (int j=0; j<n_components; j++){
      for (int k=0; k<n_components;k++){
        XHtcurr_summ+=Xmat[i][k]*Ht[k][j];
      }
      XHt[i][j]= XHtcurr_summ;
      XHtcurr_summ=0;
    }
  }
  //to multiply X * the matrix we hav eto deal with a vector and stop in time

  //printf("multiplicative matrix Xht");
  //print_matrix(XHt, X_rows,n_components);

  double violation = 0;
  double pg = 0;
  double hess = 0;
  double max_guess =0.0;
  double guess = 0.0;

  for (int i =0; i<n_components; i++){  //i in permutations
    for (int j=0; j<X_rows; j++){ //j in number of samples
      //gradient  grad = -XHt[i, t]  i in python is j here and t is i
      double grad = -XHt[j][i];

      for (int r=0; r<n_components; r++){
        //grad += HHt[t, r] * W[i, r]
        grad+=HttHt[i][r]*W[j][r];
        printf("Grad value %.2f\n", grad);
      }
      //projected gradient
      if(W[j][i]==0 && grad>0.0){
        printf("W[j][i]: %.2f ==0 and grad %.2f>0\n", W[j][i], grad);
        pg = 0;
        printf("%.2f\n",pg);
      }
      else if (W[j][i]==0 && grad<0.0){
        printf("W[j][i]: %.2f ==0 and grad %.2f<0\n", W[j][i], grad);
        pg = grad;
        printf("%.2f\n",pg);
      }
      else{
        printf("W[j][i]: %.2f !=0 and grad %.2f\n", W[j][i], grad);
        pg =grad;
        printf("%.2f\n",pg);
      }
      printf("Violation %.2f\n", violation);
      violation += fabs(pg);
      printf("Violation %.2f\n", violation);
      //Hessian
      hess = HttHt[i][i];
      //compute W if the hessian is !=0
      if (hess!=0){
        guess = (W[j][i] - grad)/hess;
        printf("Elmeents W[%d][%d] is %.2f\n",j,i, W[j][i]);
        printf("THe guess is %.2f\n", guess);
        if (guess>0.0){
          printf("The guess is >0.0 so max_guess ");
          max_guess = guess;
          printf("%.2f\n", max_guess);
        }
        else{
          printf("THe guess is <0.0 so max_Guess ");
          max_guess = 0.0;
          printf("%.2f\n", max_guess);
        }
        //update Wf
        printf("Updating matrix W[j][i] %d, %d with max_guess %.2f\n", j, i, max_guess);
        W[j][i] = max_guess;
      }
    }
  }
  printf("Updated matrix");
  print_matrix(W, X_rows, n_components);
  //print_matrix(W, X_rows, n_components);
  return violation;
}


void coordinate_descent_solver(double **H, double** W, double** Xmat, int n_components, int X_rows, int X_columns, int max_iter ){
  //in this solver   The objective function is minimized with an alternating minimization of W
  //and H. Each minimization is done with a cyclic (up to a permutation of the
  //features) Coordinate Descent.

  //use the transpose of H
  //H is n_componets * n_features, where n_features is the number of columns in X
  double violation = 0.0;
  double violation_init=0.0;

  double **Ht;
  Ht = (double **) malloc(X_columns * sizeof(double*)) ;
  for(int i = 0; i<X_columns; i++) {
        Ht[i] = (double *) malloc(n_components * sizeof(double));
  }
  double **XmatT;
  XmatT = (double **)malloc(X_columns*sizeof(double));
  for (int i= 0; i<X_columns;i++){
      XmatT[i]=(double *)malloc(X_rows *sizeof(double));
  }
  XmatT = transpose_matrix(Xmat, XmatT, X_columns, X_rows);
  //transpose
  printf("transposing matrix");
  Ht = transpose_matrix(H, Ht, X_columns, n_components);

  print_matrix(Ht, X_columns, n_components);

  for (int i=1; i< max_iter+1; i++){
      //update W
      violation+=update_coordinate_descent(Xmat,Ht, W,n_components,X_rows, X_columns);

      printf("Updated W\n");
      print_matrix(W,X_rows, n_components);
      //printf("Violation after H update %.2f\n", violation);
      //update H
      violation+=update_coordinate_descent(XmatT, W, H, n_components, X_columns, X_rows);
      printf("Updated H\n");
      print_matrix(H, n_components, X_columns);
      if (i==1){
        violation_init=violation;
      }
      if (violation_init==0){
        break;
      }
      printf("Violation/Violation_init: %.4f\n", (violation/violation_init));
      if  (violation/violation_init<tol){
        printf("Converged at iteration %d\n", (i+1));
      }
  }
  free(Ht);



}
double* nmf(int n_components, int max_iter, int random_state, double alpha, double *X,
            int X_rows, int X_columns)
{

  //transform X into a matrix
  double **Xmat;
  Xmat = (double **) malloc(X_rows * sizeof(double*));
  for (int i =0; i< X_rows; i++){
    Xmat[i] = (double *)malloc(X_columns*sizeof(double));
  }
  int counter =0 ;
  for (int i=0; i<X_rows;i++){
    for(int j=0; j<X_columns;j++){
      Xmat[i][j]=X[counter];
      counter+=1;
    }
  }
  print_matrix(Xmat, X_rows, X_columns);
  //take the average elmeent
  double avg_elem;
  avg_elem = define_avg_element(Xmat, X_rows, X_columns, n_components);
  //initialize random H
  double **H;
  H = (double **) malloc(n_components * sizeof(double*)) ;
  for(int i = 0; i<n_components; i++) {
        H[i] = (double *) malloc(X_rows * sizeof(double));
  }
  H  = initialize_H(H, avg_elem, X_rows, X_columns, n_components, 42);
  printf("Initialized H matrix");
  print_matrix(H, n_components, X_columns );
  //initialize random W
  double **W;
  W = (double **)malloc(X_rows *sizeof(double*));
  for (int i=0; i<X_rows; i++){
      //X_rows = n_samples   X_columns = n_features
      W[i]=(double *) malloc(n_components *sizeof(double));
  }
  W = initialize_W(W, avg_elem, X_rows, X_columns, n_components, 42);
  printf("Initialized W matrix");
  print_matrix(W, X_rows, n_components);
  printf("Coordinate descent\n");
  //TODO regularization, skip for the moment
  //Coordinate descent solver
  coordinate_descent_solver(H, W, Xmat, n_components, X_rows, X_columns,max_iter);

  return H, W;
  free(H);
  free(W);
}
