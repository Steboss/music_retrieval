
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <time.h>
#include "polyfit.h"
//#include "dfa.h"

#define M_PI 3.14159265358979323846
#define PI 3.14159265
#define LEN(arr) ((sizeof (arr) / sizeof (*arr))) //macro to define the size of an array
#define tol 0.0001


double Log2( float n )
{
    // log(n)/log(2) is log2.
    return log( n ) / log( 2 );
}

int floor_div(int a, int b)
{
    assert(b != 0);
    div_t r = div(a, b);
    if (r.rem != 0 && ((a < 0) ^ (b < 0)))
        r.quot--;
    return r.quot;
}
/**************/
float vector_avg(float *X, int n_elems){
  float avg = 0;
  for (int i=0; i<n_elems;i++){
    avg+=X[i];
  }
  avg /= n_elems;
  return avg;
}


float* cumsum_averaged(float* X, int n_elems, float avg, float* Xcumsum){
  /*Cumulative sum function. Sum each element of a vector and return a cumsum
  vector of the same dimensions as the input X (n_elems)
  E.g.
  X
  [1,2,3,4,5,6]
  cumsum:
  [1,3,6,10,15,21]

  Parameters
  ---------
  X: float *
     input X signal
  n_elems: int
     number of elements of X
  avg: float
     average of signal X
  Xcumsum: float *
    cumulative sum vector of X

  */

  for (int i=0; i<n_elems;i++){
    if (i==0){
      Xcumsum[i] = X[i]-avg;
    }
    else{
      Xcumsum[i] = Xcumsum[i-1] +(X[i]-avg);

    }
  }

  return Xcumsum;
}



void dfa(float* X, int n_elems, int* scales, int scales_len, float* fluct, float* coeffs){
  /*dfa function  -- at the moment let's stick to a vectorial view, then we can think
  of having amulti dimensional dfa with dataframes/array
  Parameters
  ----------
  X: float*
     input signal vector
  n_elems: int
     number of elements in vector X
  scales: float*
     vector of scales, start from scale_low to scale_max with scale_dens as range
    created in python -- explanation here:
                      scale_low: int
                         first scale value
                      scale_max: int
                         last scale value
                         scale values will be converted to an interval range
                         from 2^scale_low up to 2^scale_max
                      scale_dens; float
                         how dense should be the scale_low -- scale_max interval
                         range. E.g.:
                         scale_low = 5, scale_max = 9, scale_dens=0.25 :
                         intervalsfrom 2^5 to 2^9 with 0.25
                         step between values

  Return
  ------

  */
  //fix the fitting order
  int fitting_order = 1;
  float* Xcumsum;
  Xcumsum = (float*)malloc(n_elems*sizeof(float));
  //average of the vector X
  float avg = vector_avg(X, n_elems);
  //then cumulative sum - average  --> y = np.cumsum(x - np.mean(x))
  Xcumsum = cumsum_averaged(X, n_elems, avg, Xcumsum);

  //final results will be stored in rms vector
  //rms vector has the dimension of len of scale --> scales_len
  float * rms; //vector whose dimensions are  as the same as xcumsum
  rms = (float*) malloc(scales_len*sizeof(int));
  //initialize to 0
  for (int i=0; i<scales_len;i++){
      rms[i]= 0.0;
  }

  //initialize the coefficient vector
  //initialize a vector whose shape is fitting_Order +1 to store the coefficient of the fitting
  float* tmpCoeff;
  tmpCoeff = (float*)malloc((fitting_order+1)*sizeof(float)); //2 is 1 + 1 where 1 is the order of th epolynomiaal fitting
  for (int h=0; h<(fitting_order+1);h++){
    tmpCoeff[h]=0.0;
  }//if we had order 3, we would have had 4 as malloc(4*sizeof(float))

  //HERE WE START THE RMS calculation
  for (int i=0; i<scales_len; i++){
      //take the current scale
      int curr_scale = scales[i];
      printf("Current scale %d\n", curr_scale);
      //create a scale_ax vector
      //compute the shape
      int shape1 = floor_div(n_elems, curr_scale);
      int counter =0  ; //this cound will be updated in order to cycle thorugh all the elmeents
                      //of x in correct order, updating the strided single vecto relment tmpX
      //initialize a tmporary X vector, to host single strided Xcumsum matrix
      float* tmpX;
      tmpX = (float*) malloc( curr_scale*sizeof(float)); //this will contain a max of curr_scale elements
      //initialize to 0
      for (int h=0; h<curr_scale;h++){
        tmpX[h] = 0.0;
      }
      //now we have the first bit of X with curr_scale elements
      //we have to perform the fitting so we need an x-axis, which is made of
      //number of elements of the tmpX vectors, so curr_scale number of elemetns

      //basically it goes from 0 to curr_scale - 1
      //now polyfit of scale_ax and tmp  --> extenral funciton
      float* scale_ax;
      scale_ax = (float*)malloc(curr_scale*sizeof(float*));
      //
      for(int h=0; h< curr_scale;h++){
        scale_ax[h] = h;
      }
      //
      //initialize a vector of size curr_scale, as tmpX, to be used to evaluate
      //the goodness of fitting  between fitting_result and tmpX[counter2] element
      float* tmpDiff;
      tmpDiff = (float*)malloc((curr_scale)*sizeof(float));
      for (int h=0; h<(curr_scale);h++){
        tmpDiff[h]=0.0;
      }
      //
      //store all the rms in a shape1 vector
      float* tmpShape1;
      tmpShape1=(float*)malloc( (shape1)*sizeof(float));
      for (int h=0; h<shape1;h++){
        tmpShape1[h]=0.0;
      }

      //now we have to divide Xcumsum into shape1 elements, each of which has curr_scale elements
      //whcih are updated, namely move along the Xcumsum vector, thorugh counter
      for (int j=0; j< shape1; j++){ //this is the number of final matrix we want to have
          for (int k=0; k< curr_scale;k++){ //this is the number of the elements of X to be processed at each j step
                tmpX[k] = Xcumsum[counter];
                counter+=1;
            }
          //POLYFIT
          //perform the fitting and update tmpCoeff
          polyfit(scale_ax, //int*
                  tmpX, //float*
                  curr_scale, //int
                  fitting_order, //try a linear fitting
                  tmpCoeff //float*
                ) ;
          printf("tmpCoeff[0] :%.4f  tmpCoeff[1] : %.4f\n", tmpCoeff[0], tmpCoeff[1]);
          //now we need to evaluate the polynomial fit
          //to do that we can do the  polyval so  np.polyval(coeff, scale_ax)
          //if we have two coefficient we have  tmpCoeff[1]*x + tmpCoeff[0]
          //x must be substituted with the i-th element of tmpX
          float fitting_result = 0.0;
          //
          for(int t=0; t<curr_scale;t++){
              for(int u=(fitting_order); u>-1;u--){
                  if(u==0){
                    //this is the last step
                    fitting_result+= tmpCoeff[u]; //here no tmpX[t] because we have 0 for x order so x^0=1
                  }
                  else{
                    //tmpCoeff * x^fitting_order - u
                    //e.g. tmpCoeff[0]*x  --> fitting-Order = 1 and u = 0
                    fitting_result+=tmpCoeff[u]*(pow(scale_ax[t],u));
                  }
              }
              //now here we can compute how good is the fit between element t and element
              //t of the tmpX at the moment we are computing the difference of elements
              float curr_diff = tmpX[t]-fitting_result;
              float squared_diff =  pow( curr_diff,2);
              printf("tmpX[%d]: %.4f, fitting_result:%.4f, diff:%.4f, pow:%.4f\n",
                     t, tmpX[t], fitting_result, curr_diff, squared_diff);
              tmpDiff[t] =squared_diff;
              printf("pow(tmpX[t] - fitting_result, 2):%.8f\n", tmpDiff[t]);
              fitting_result = 0.0;
              //at the end we need to compute average of ((xcut-xfit)**2))
          }
          //here we are closing shape1 cycle, namely cycle through all the strided elements of X
          //at this stage we need to compute the sqrt of the average of tmpShape1 elements
          float tmpDiffavg = vector_avg(tmpDiff, curr_scale);
          //this is the shape1-th element of the final
          tmpShape1[j] = sqrt(tmpDiffavg);
          printf("elment %d  rms %.8f\n", j, tmpShape1[j]);
      }
      //
      //printf("Freeing tmpDiff\n");
      free (tmpDiff);
      tmpDiff = NULL;
      //
      ///printf("Freeing tmpX\n");
      free (tmpX);
      tmpX = NULL;
      //
      //printf("Computing sqrt of average vector of tmpShape1\n");
      for (int h=0; h<shape1;h++){
        tmpShape1[h] = pow(tmpShape1[h],2);
      }
      float tmpShape1avg = vector_avg(tmpShape1, shape1);
      tmpShape1avg = sqrt(tmpShape1avg);
      //printf("Final result %.4f\n", tmpShape1avg);
      //this is the rms i-th sotion s
      //Check values in rms
      //printf("Updating element i %d oftmprms vector\n",i);
      rms[i] = tmpShape1avg;
      printf("Fluctuation %d : %.4f\n", i, rms[i]);
      //printf("Freeing tmpShape1\n");
      free (tmpShape1);
      tmpShape1 = NULL;
  }

  //What we have computed are the fluctuations
  //now we can retrieve the Hurst exponent via a polyfit again
  float* log2_scales ;
  log2_scales = (float*)malloc(scales_len*sizeof(float));
  for(int h=0; h<scales_len;h++){
    log2_scales[h]= Log2( (float)scales[h]);
  }


  float* log2_fluct ;
  log2_fluct = (float*)malloc(scales_len*sizeof(float));
  for(int h=0; h<scales_len;h++){
    log2_fluct[h]= Log2(rms[h]);
  }

  //1 coefficient
  polyfit(log2_scales,
          log2_fluct,
          scales_len,
          fitting_order,
          coeffs
        );

  for (int h=0; h<scales_len;h++){
    fluct[h] = rms[h];
  }
  //after this we'll need to compute the fitting of the rms
  //printf("Freeing tmpCoeff\n");
  free (tmpCoeff);
  //printf("Freeing cumsum\n");
  free(Xcumsum);
  //printf("Freeing rms\n");
  free(rms);
  //
  free(log2_fluct);
  free(log2_scales);
}
