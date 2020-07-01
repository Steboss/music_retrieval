//One thing to be checked is if stft needs to be raised at te power of 2
//or cabs is enough. Power = 1 .. energy, Power= 2 Power

/*Idea:
1. Compute hte stft of the input signal  and retrieve the magnitude
2. Create a mel filter --> use the stft frequencies, translate to mel_frequencies --> compute weights
3. dot product between stft magnitude and mel_basis (2.)
4. compute the discrete cosine fourier transform with fftw andd retrieve results as [:n_mfcc]
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <time.h>

#include "stft.h"

#define M_PI 3.14159265358979323846
#define PI 3.14159265
#define LEN(arr) ((sizeof (arr) / sizeof (*arr))) //macro to define the size of an array
#define tol 0.000001
//HELP FUNCTIONS
double* linspace(double *x, double init, double fin, int N){
  //create a linspace array
  double step = (fin - init)/(N);
  x[0] = init;
  for (int i=0; i<N; i++){
    x[i] = x[i-1]+step;
  }
  x[N-1] = fin;
  return x;
}

////////////////////////////////////////////////////////////////////////////
double mel_to_hz(double curr_mel){
  //convert mel frequency to hz
  double f_min = 0.0 ;
  double f_sp = 200.0/3.0;
  // nonlinear scale
  double min_log_hz =1000.0;  //beginning of log region (Hz)
  double min_log_mel = (min_log_hz - f_min)/f_sp;  //same (mels)
  double logstep = log(6.4)/27.0 ; //step size for log region
  //double freqs = 0.0;
  if (curr_mel >=min_log_mel){
    return  min_log_hz * exp(logstep*(curr_mel-min_log_mel));
  }
  else{
    return  f_min + (f_sp*curr_mel);
  }
}

double hz_to_mel(double frequency){
  //convert the frequency to mel scale
  double f_min = 0.0;
  double f_sp = 200.0/3.0 ;
  //log scale part
  double min_log_hz = 1000.0 ;
  double min_log_mel =(min_log_hz - f_min)/f_sp;
  double logstep = log(6.4)/27.0 ;
  //double mels = 0.0;
  if (frequency>= min_log_hz){
    return   min_log_mel + log(frequency/min_log_hz)/logstep;
  }
  else{
    return   (frequency - f_min)/f_sp ;
  }
}
double* mel_frequencies(int n_mels, double fqmin, double fqmax, double* melfreqs){
  /*acoustic frequencies tuned to the mel scale.
    The mel scale is a quasi-logarithmic function of acoustic frequency
    designed such that perceptually similar pitch intervals (e.g. octaves)
    appear equal in width over the full hearing range.
    Because the definition of the mel scale is conditioned by a finite number
    of subjective psychoaoustical experiments, several implementations coexist
    in the audio signal processing literature [1]_. By default, librosa replicates
    the behavior of the well-established MATLAB Auditory Toolbox of Slaney [2]_.
    According to this default implementation,  the conversion from Hertz to mel is
    linear below 1 kHz and logarithmic above 1 kHz. Another available implementation
    replicates the Hidden Markov Toolkit [3]_ (HTK) according to the following formula:
    `mel = 2595.0 * np.log10(1.0 + f / 700.0).`
 */
    double min_mel = hz_to_mel(fqmin);
    double max_mel = hz_to_mel(fqmax);
    linspace(melfreqs, min_mel, max_mel, n_mels);
    //convert back to hz
    for (int i=0; i<n_mels; i++){
      melfreqs[i] = mel_to_hz(melfreqs[i]);
    }

    return melfreqs;
}

void mel(int sample_freq,int windowSize, double** weights ){
  //Create a Filterbank matrix to combine FFT bins into Mel-frequency bins
  //TODO one parameter to check is n_mels. Here we hard coded to 128
  //TODO fmin and fmax are hard coded as well
  // use Slaney formula
  double fqmin = 0.0 ;
  double fqmax = (sample_freq/2) ;//fqmax can be user-defined, so you can decide to crop up to fqmax
  int n_mels = 128;
  //create a matrix for weights
  int n_cols = 1+floor(windowSize/2);

  printf("n_cols %d\n",n_cols);
  //center frequencies
  printf("Defining fftfreqs\n");
  double *fftfreqs;
  fftfreqs = (double*)malloc(sizeof(double)*(n_cols));  //n_cols vector
  linspace(fftfreqs, fqmin, fqmax, n_cols);
  //center frequencies of mel bands uniformly
  printf("Defining melfreqs\n");
  double *melfreqs;
  int melfreqs_size = n_mels+2;
  melfreqs = (double*)malloc(sizeof(double)*(melfreqs_size));
  mel_frequencies(melfreqs_size, fqmin, fqmax, melfreqs);

  //nmels_+2 vector

  //compute the n-th discrete differnece in the mel frequencies
  //first diffence  melfreqs[i] = melfreqs[i+1] - melfreqs[i]
  //and so on
  printf("Defining difference\n");
  double *fdiff;
  fdiff = (double*)malloc(sizeof(double)*(n_mels+1)); //one element less than mefreqs
  for(int i=0; i<n_mels+1; i++){
    fdiff[i] = melfreqs[i+1] - melfreqs[i];
  }
  //n_mels+1 vector
  //outer subtraction
  //the result will be a matrix whose dimensions are n_cols+2 * n_cols
  //each element will be the subtraction of a fixed mel_f element against all the fftfreqs
  /* np.subtract.outer([1,2,3],[4,5])
  Out[23]:
  array([[-3, -4],
        [-2, -3],
        [-1, -2]])
  */
  printf("Defining ramps\n");
  int n_mels_2 = n_mels + 2;
  double **ramps ; //dimensions n_mels+2 * n_cols
  ramps = (double**) malloc(n_mels_2*sizeof(double*));
  for (int i=0; i<(n_mels_2);i++){
    ramps[i] = (double*)malloc(n_cols*sizeof(double));
  }

  for (int i=0; i<n_mels+2;i++){
    for (int j=0; j<n_cols; j++){
      ramps[i][j] = melfreqs[i] - fftfreqs[j];
    }
  }

  //lower and upper slopes for all bins
  //double* lower_array ;
  //lower_array = (double*)malloc(sizeof(double)*(n_cols));//vector n_cols dimension
  //double * upper_array ;
  //upper_array = (double*)malloc(sizeof(double)*(n_cols)); //vector n_cols dimensions
  double lower_val = 0.0;
  double upper_val = 0.0;
  double minimum_val = 0.0;
  //define row eleements for weights by column
  //now divide all the ramps by the elements in fdiff
  for(int i=0; i<n_mels; i++){
    for (int j=0; j<n_cols;j++){
      lower_val = (-ramps[i][j])/fdiff[i];
      upper_val = (ramps[i+2][j]/fdiff[i+1]);//is this right!?
      //construct here the minumum array
      if (lower_val>upper_val){
        minimum_val = upper_val;
      }
      else{
        minimum_val = lower_val;
      }
      //now fill the weight
      if (minimum_val>0.0){
        weights[i][j] = minimum_val;
      }
      else{
        weights[i][j] = 0.0;
      }
      //printf("weights %.2f\n", weights[i][j]);
    }
  }
  printf("Defining enorm");
  //slaney norm
  double *enorm;
  enorm = (double*)malloc(sizeof(double)*(n_mels)); //I think here we're retrieving the 2 last elems
  for (int i=0; i<n_mels;i++){
    enorm[i] = 2.0/(melfreqs[i+2]- melfreqs[i]);
  }
  //multiply weight with enorm
  for (int i=0; i<n_mels;i++){
    for (int j=0; j<n_cols;j++){
      weights[i][j]*=enorm[i] ; //TODO DOUBLE CHECK THIS CALCUALTION!
    }
  }

  //return weights;

  free(fftfreqs);
  free(melfreqs);
  free(fdiff);
  for (int i=0; i< n_mels+2; i++){
    free(ramps[i]);
  }
  free(ramps);
  free(enorm);


}

double* mfcc(double *wav_data, int samples, int windowSize, int hop_size,\
             double *magnitude, int sample_freq, int length, int n_mfcc, double *cython_basis)
{
    //1. Compute stft of the input signal
    //we could pass this info from cython already as a matrix
    //stft(wav_data, samples, windowSize, hop_size, magnitude, sample_freq, length);


    double** magnitude_array ;
    int cols = ( (length/(windowSize/2)) -1);  //109
    int rows = 1+floor(windowSize/2);  //1025
    printf("C dimensions: %d, %d, %d\n", cols, rows, cols*rows);
    magnitude_array = (double**)malloc(rows*sizeof(double*));
    for(int i=0; i<rows;i++){
      magnitude_array[i] =(double*)malloc(cols*sizeof(double));
    }

    int counter = 0 ;
    for (int i =0; i< rows; i++){
      for (int j=0; j< cols; j++){
        magnitude_array[i][j] = magnitude[counter];
        counter+=1;
      }
    }
    printf("Magnitude transposed to array\n");

    //2. mel filter
    int n_mels = 128; //TODO this is hard coded  //128
    int n_cols = 1+floor(windowSize/2);  //1025
    double **weights ;
    weights = (double**) malloc(n_mels*sizeof(double*));
    for (int i=0; i<n_mels;i++){
      weights[i] = (double*)malloc(n_cols*sizeof(double));
    }

    for (int i=0; i<n_mels;i++){
      for(int j=0; j<n_cols;j++){
        weights[i][j]=0.0;
      }
    }

    printf("Mel filter\n");
    mel(sample_freq,windowSize, weights ) ; //the output matrix should be
    // M         : np.ndarray [shape=(n_mels, 1 + n_fft/2)]
    //3.dot product
    /*STFT magnitude  windowsSize/2 + 1  x  length(windowSize*0.5) -1
    / M               n_mels, 1+n_fft/2
    //Thus the product M * STFT  can be done as we have
    // n_mels * 1+n_fft/2   --- n_fft/2 + 1 x length...*/

    //double check the weights

    printf("\nDot product n_mels * stft\n");

    //M  n_mels * 1+n_fft/2
    //S 1+n_fft/2 *  length*(n_fft*0.5)-1
    //result   n_mels * lengt*(n_fft*0.5)-1
    //initialize the  mel_basis result

    //magnitude_array   (1+ nnft/2) * (length....)
    //M                 n_mels * 1+nnft/2
    //M* magnitude_array  n_mels*length
    int basis_cols =  ( (length/(windowSize/2)) -1);  //109
    int basis_rows = n_mels;  //128

    double **mel_basis = (double**)malloc(basis_rows*sizeof(double*));
    for (int i=0 ; i<basis_rows;i++){
      mel_basis[i] = (double*)malloc(basis_cols*sizeof(double));
    }
    printf("initialize mel basis\n");
    for (int i=0; i< basis_rows;i++){
      for (int j=0; j<basis_cols;j++){
        mel_basis[i][j] = 0.0;
      }
    }
    /*remember matrix moltiplication:
    int a[M][N]; == weights   n_mels *  nnft
    int b[N][P];==  magnitude nnft *  length
    int c[M][P]; == mel_basis n_mels * length
    MPN  --> for cycle, n_mels, length, nfft*/
    double log_spectro = 0.0;
    double amin = 0.00000001;
    //double ref_value = 1.0; //this is 1 so the log will be = 0.0
    double top_db = 80.0; //max convertible decibel
    double log_spectro_max = 0.0; //
    printf("Compute mel_basis product\n");
    for (int i=0; i< n_mels; i++){ //128
      for (int j=0; j< cols; j++){  //109
        for (int k=0;k< rows; k++){  //128
          //printf("weights[%d][%d]:%.2f\n",i, k, weights[i][k]);
          //printf("magnitude_array[%d][%d]:%.2f\n", k,j, magnitude_array[k][j]);
          mel_basis[i][j] += weights[i][k]*magnitude_array[k][j];
          //printf("i:%d/%d, j:%d/%d, k:%d/%d\n", i, n_mels, j, cols, k,rows);
          //printf("weights[i][j]: %.2f, magnitude_array[k][j]: %.2f\n", weights[i][k], magnitude_array[k][j]);
        }
        //convert the element to dB from amplitude square
        //convert to dB
        if (mel_basis[i][j]>amin){
          log_spectro = 10.0*log10(mel_basis[i][j]);
          log_spectro_max = log_spectro - top_db;
        }
        else{
          log_spectro = 10.0*log10(amin);
          log_spectro_max  = log_spectro - top_db;
        }
        //printf("log_spectro: %.2f, log_spectro_max: %.2f\n", log_spectro, log_spectro_max);
        if (log_spectro> log_spectro_max){
          mel_basis[i][j] =log_spectro;
          //printf("log_spectro> log_spectro_max");
        }
        else{
          mel_basis[i][j] =log_spectro_max;
          //printf("log_spectro< log_spectro_max");
        }
        //printf("mel_basis[%d][%d]:%.2f\n",i,j,mel_basis[i][j]);

      }
    }

     //4. COmpute the discrete cosine fourier transform and take all the
     //element till n_mfcc, namely the number of parameters we want to extract
    counter = 0 ;
    for (int i=0; i<basis_rows;i++){
      for (int j=0; j<basis_cols;j++){
        cython_basis[counter] = mel_basis[i][j];
        counter+=1;
      }
    }
     //fftw_plan plan = fftw_plan_r2r_1d(X,X,X, FFTW_REDFT10, FFTW_ESTIMATE);
     //fftw_execute(plan);
    for (int i=0; i< rows; i++){
         free(magnitude_array[i]);
    }
    free(magnitude_array);
    for (int i=0; i< n_mels; i++){
         free(weights[i]);
    }
    free(weights);
    for (int i=0; i<basis_rows;i++){
     free(mel_basis[i]);
    }
    free(mel_basis);

    return cython_basis;


}
