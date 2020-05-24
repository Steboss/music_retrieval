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

float ** multiply_matrix(float **mat, float **mat2, float ** result, int row, int col){

  float ** backup = result;
  for (int i =0; i< row; i++){
    for(int j =0; j< col; j++){
      for (int k=0; k< col; k++){
        backup[i][j]+= mat[i][k]*mat2[k][j];
      }
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
      summa+=X[i][j];
    }
  }
  if (summa<0.0){
    summa = - summa;
  }
  float avg = sqrt( (summa/n_elems)/ n_components) ;

  printf("%.2f\n",avg);
  return avg;
}

float** initialize_H(float** H, float avg_elem, int X_rows, int X_columns, int n_components, int r_seed){
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

float** initialize_W(float ** W,
                      float avg_elem,
                      int X_rows,
                      int X_columns,
                      int n_components,
                      int r_seed){

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

float** coordinate_descent_solver(float **H,
                                   float** W,
                                   float** Xmat,
                                   int n_components,
                                   int X_rows,
                                   int X_columns,
                                   int max_iter,
                                   int verbose ){
  //in this solver   The objective function is minimized with an alternating minimization of W
  //and H. Each minimization is done with a cyclic (up to a permutation of the
  //features) Coordinate Descent.

  //use the transpose of H
  //H is n_componets * n_features, where n_features is the number of columns in X
  float violation = 0.0;
  float violation_init=0.0;
  float violation_middle=0.0;
  float pg = 0;
  float hess = 0;
  float max_guess =0.0;
  float guess = 0.0;
  float grad = 0.0;
  //Transpose of H
  if (verbose==1){
    printf("Initializing matrix for coordinate descent\n");
  }
  float **Ht;
  Ht = (float **) malloc(X_columns * sizeof(float*)) ;
  for(int i = 0; i<X_columns; i++) {
        Ht[i] = (float *) malloc(n_components * sizeof(float));
  }
  Ht = transpose_matrix(H, Ht, X_columns, n_components);
  if (verbose==1){
    printf("H transposed!\n");
  }
  //Transpose of W
  float **Wt;
  Wt = (float **) malloc(n_components * sizeof(float*)) ;
  for(int i = 0; i<n_components; i++) {
        Wt[i] = (float *) malloc(X_rows * sizeof(float));
  }
  Wt = transpose_matrix(W, Wt, n_components, X_rows);
  if (verbose==1){
    printf("W transposed!\n");
  }
  //Transpose of X
  float **XmatT;
  XmatT = (float **)malloc(X_columns*sizeof(float*));
  for(int i = 0; i<X_columns; i++) {
        XmatT[i] = (float *) malloc(X_rows * sizeof(float));
  }
  XmatT = transpose_matrix(Xmat, XmatT, X_columns, X_rows);
  if (verbose==1){
    printf("X transposed!\n");
  }
  //Get an empty 2 D array
  float **HtT;
  HtT = (float ** )malloc(n_components*sizeof(float*));
  for (int i=0; i<n_components; i++){
    HtT[i] = (float*)malloc(X_columns*sizeof(float));
  }
  if (verbose==1){
    printf("HtT done!\n");
  }
  //
  float **HttHt;
  HttHt = (float**) malloc(n_components*sizeof(float*));
  for (int i=0;i<n_components;i++){
    HttHt[i]=(float*)malloc(n_components*sizeof(float));
  }
  if (verbose==1){
    printf("HttHt done!\n");
  }
  //
  float **XHt;
  XHt = (float**) malloc(X_rows*sizeof(float*));
  for (int i=0;i<X_rows;i++){
    XHt[i]=(float*)malloc(n_components*sizeof(float));
  }
  if (verbose==1){
    printf("XHt done!\n");
  }

  //For updating H
  float **WttWt;
  WttWt = (float**) malloc(n_components*sizeof(float*));
  for (int i=0;i<n_components;i++){
    WttWt[i]=(float*)malloc(n_components*sizeof(float));
  }//
  if (verbose==1){
    printf("WttWt done!\n");
  }
  //
  float **XtW;
  XtW = (float**) malloc(X_columns*sizeof(float*));//
  for (int i=0;i<X_columns;i++){ //X_columns
    XtW[i]=(float*)malloc(n_components*sizeof(float));
  }
  if (verbose==1){
    printf("XtW done!\n");
  }

  if (verbose==1){
    printf("Coordinate descent initialized\n");
  }

  ///iterationsf
  for (int i=1; i< max_iter+1; i++){
      if (verbose==1){
        printf("Iteration: %d\n", i);
      }
      violation=0.0;
      //update W
      //=update_coordinate_descent(Xmat,Ht, W,n_components,X_rows, X_columns);
      //1. transpose Ht
      //printf("updating W\n");
      HtT = transpose_matrix(Ht, HtT, n_components, X_columns);
      //printf("HtT matrix\n");
      //print_matrix(HtT, n_components, X_columns);
      //2. Dot product HtT and HT
      float HttHtcurr_summ = 0.0;
      for (int i=0; i<n_components; i++){
        for(int j=0; j<n_components;j++){
          for(int k=0; k< X_columns;k++){
            HttHtcurr_summ+=HtT[i][k]*Ht[k][j];
          }
          //fill the array
          HttHt[i][j] = HttHtcurr_summ;
          HttHtcurr_summ=0;
        }
      }
      //printf("HttHt matrix\n");
      //print_matrix(HttHt, n_components, n_components);
      //3. dot X and Ht --> X-rows *n_components
      float XHtcurr_summ= 0.0;
      for (int i =0; i<X_rows; i++){
        for (int j=0; j<n_components; j++){
          for (int k=0; k<X_columns;k++){
            XHtcurr_summ+=Xmat[i][k]*Ht[k][j];
          }
          XHt[i][j]= XHtcurr_summ;
          XHtcurr_summ=0;
        }
      }
      //printf("XHT matrix\n");
      //print_matrix(XHt, X_rows, n_components);
      //4.permutations to update W
      violation_middle = 0.0;
      //printf("W computation");
      for (int i =0; i<n_components; i++){  //i in permutations --> t

        for (int j=0; j<X_rows; j++){ //j in number of samples  --> i
          //gradient  grad = -XHt[i, t]  i in python is j here and t is i
          grad = -XHt[j][i];
          for (int r=0; r<n_components; r++){  //--> r
            //grad += HHt[t, r] * W[i, r]
            grad+=HttHt[i][r]*W[j][r];
          }
          //projected gradient
          if(W[j][i]==0 && grad>0.0){
            pg = 0;
          }
          else if (W[j][i]==0 && grad<0.0){
            pg = grad;
          }
          else{
            pg =grad;
          }
          violation_middle += fabs(pg);
          //Hessian
          hess = HttHt[i][i];

          //compute W if the hessian is !=0
          if (hess!=0){
              guess = (W[j][i]) - (grad/hess);
              if (guess>0.0){
                max_guess = guess;
              }
              else{
                max_guess = 0.0;
              }
              //update Wf
              W[j][i] = max_guess;
          }

        }
      }
      violation +=violation_middle;
      //printf("Updated  W matrix:\n");
      //printf("UPDATE W MATRIX");
      //print_matrix(W, X_rows, n_components);
      //update W transpose!!!!
      Wt = transpose_matrix(W, Wt, n_components, X_rows);
      //Update H
      //1. Transpose of W --> WT -- Wt
      //2. dot product WT and W  [2,6]*[6,2]
      //printf("Updating H\n");
      float WttWtcurr_summ = 0.0;
      //W has dimensions 6,2  and Wt  2, 6
      //we are moltipliiyng Wt X W  so we'll get a n_components * n_compoments matrix
      for (int i=0; i<n_components; i++){
        for(int j=0; j<n_components;j++){
          for(int k=0; k< X_rows;k++){
            WttWtcurr_summ+=Wt[i][k]*W[k][j];
          }
          //fill the array
          WttWt[i][j] = WttWtcurr_summ;
          WttWtcurr_summ=0;
        }
      }
      //printf("WttWt matrix");
      //print_matrix(WttWt, n_components, n_components);
      //3. dot Xt and W
      float XWtcurr_summ= 0.0;
      for (int i =0; i<X_columns; i++){//X_Columns
        for (int j=0; j<n_components; j++){
          for (int k=0; k<X_rows;k++){ //n_components
            XWtcurr_summ+=XmatT[i][k]*W[k][j];
          }
          XtW[i][j]= XWtcurr_summ;
          XWtcurr_summ=0;
        }
      }
      //printf("XWt matrix");
      //print_matrix(XtW, X_columns, n_components);
      //1025 x 2  or  95 x 2
      //4. permutations
      //original W, HHT, XHT
      //Here we have Ht, WttWt, XWt
      //n_components here derived from Ht dimensions
      //the second look is through the n_samples which is Ht.sahep[0
      violation_middle =0.0;
      for (int i =0; i<n_components; i++){  //i in permutations
        for (int j=0; j<X_columns; j++){ //j in number of samples
          //gradient  grad = -XHt[i, t]  i in python is j here and t is i
          grad = -XtW[j][i];
          for (int r=0; r<n_components; r++){
            grad+=WttWt[i][r]*Ht[j][r];
          }
          //projected gradient
          if(Ht[j][i]==0 && grad>0.0){
            pg = 0;
          }
          else if (Ht[j][i]==0 && grad<0.0){
            pg = grad;
          }
          else{
            pg =grad;
          }
          violation_middle += fabs(pg);
          //Hessian
          hess = WttWt[i][i];
          //compute W if the hessian is !=0
          if (hess!=0){
            guess = Ht[j][i] - (grad/hess);
            if (guess>0.0){
              max_guess = guess;
            }
            else{
              max_guess = 0.0;
            }
            //update Wf
            Ht[j][i] = max_guess;
          }
        }
      }
      //print_matrix(Ht, n_components, X_columns);
      violation+=violation_middle;
      //printf("Update H matrix");
      //print_matrix(Ht, n_components, X_columns);
      //printf("Current violation %.6f\n", violation);

      if (i==1){
        violation_init=violation;
        //printf("Violation init: %.6f\n", violation_init);
      }
      if (violation_init==0){
        break;
      }
      if (verbose==1){
        printf("Violation/Violation_init: %.6f\n", (violation/violation_init));
        //printf("Current violation %.6f\n", violation);
        //printf("Current violation_init %.6f\n", violation_init);
      }
      if  (violation/violation_init<0.00001){
        printf("Converged at iteration %d\n", (i+1));
        break;
      }
  }
  //transpose H
  H = transpose_matrix(Ht, H, n_components, X_columns);

  //free memory
  free(Ht);
  free(Wt);
  free(XmatT);
  free(HttHt);
  free(XHt);
  free(WttWt);
  free(XtW);
  return 0;

}

//TODO: add the regulization variables alpha
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
