#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>

#define PI 3.14159265

//typedef double fftw_complex[2];

double hamming(int windowLength, double *buffer){
  for(int i=0; i<windowLength; i++){
    buffer[i] = 0.54 - 0.46*cos((2*PI*i)/(windowLength-1.0));
  }
  return *buffer;

}

//the main function should compute the  stft of the signal, return the stft
//and the correlation points between a chosen window --for the moment --
//and all the other windows
double* stft(double *wav_data, int samples, int windowSize, int hop_size,\
             double *magnitude, int sample_freq, int length)
{
  //samples is now the  wav file length
  printf("Initialization of parameters...\n");
  //integer indexes
  int i,counter ;
  counter = 0 ;
  //window to be applied to data
  double hamming_result[windowSize];
  double summa;
  //define the fourier variables
  fftw_complex *stft_data, *fft_result, *storage;
  //fftw_complex *fft_result, *storage;
  //makdouble *stft_data;
  //allocate the memory for the fftw
  stft_data = (fftw_complex*)(fftw_malloc(sizeof(fftw_complex)*(windowSize)));
  fft_result= (fftw_complex*)(fftw_malloc(sizeof(fftw_complex)*(windowSize)));
  storage = (fftw_complex*)(fftw_malloc(sizeof(fftw_complex)*(samples)));
  printf("Total length of storage %d\n", (samples));
  //exit(0);
  //define the fftw plane
  fftw_plan plan_forward;
  //plan_forward = fftw_plan_dft_r2c_1d(windowSize, stft_data, fft_result,FFTW_ESTIMATE);//, FFTW_FORWARD, FFTW_ESTIMATE);
  plan_forward = fftw_plan_dft_1d(windowSize,stft_data,fft_result, FFTW_FORWARD,FFTW_ESTIMATE);
  //double for the correlation value
  printf("Creation of a hamming window...");
  hamming(windowSize, hamming_result);
  //compute the scaling factors
  for (i=0; i<windowSize; i++)
  {
    summa+=hamming_result[i]*hamming_result[i];
    //printf("%.2f\n", hamming_result[i]);
  }
  //double scale = sqrt(1.0/(summa*summa));
  //now here we need to implement the stft
  int chunkPosition = 0; //actual chunk position
  int readIndex ; //read the index of the wav file
  int n_elem_read = 0 ; //number of elements read -- for sanity check only

  while (counter < samples- windowSize ){
    //read the window
    for(i=0; i<windowSize; i++){

      readIndex = chunkPosition + i;
      //printf("Index position %d\n", readIndex);
      //printf("Hamming result %d\n", hamming_result[i]) ;
      stft_data[i] = wav_data[readIndex]*hamming_result[i];//*_Complex_I  + 0.0*I;
      n_elem_read+=1;
      //printf("%.2f,",wav_data[readIndex]);

    }
    //printf("\n");
    //compute the fft
    fftw_execute(plan_forward);
    //store the half of the FFT in a data structure
    for (i=0; i<windowSize/2 +1 ;i++)
    {
      //printf("RE: %.8f  IM: %.8f\n", creal(fft_result[i]/(windowSize)),cimag(fft_result[i]/(windowSize)));
      storage[counter] = fft_result[i]; //fft_result[i];///windowSize;//creal(fft_result[i]) + cimag(fft_result[i]);
      //correlation_result[counter] = fft_result[i]*conj(fft_result[i])*scale;
      //printf("%.8f  +  %.8f\n", fft_result[i]);
      counter+=1;

    }

    chunkPosition += hop_size/2;
    //printf("Chunk Position %d\n", chunkPosition);
    //printf("Counter position %d\n", counter);
    //printf("Numb of read elements %d\n", n_elem_read);
    //printf("Numb of samples %d\n",samples);
    //printf("Length of the wav %d\n", length);
    //printf("ReadIndex %d\n", readIndex);
    //printf("Fourier transform done\n");
    //printf("-----------------------------------------------------\n");


  }
  printf("%d\n", counter);

  //printf("This is the storage capacity %d and this the counter %d\n", storage_capacity,counter);
  //normalize the fourier coefficients
  for (i=0; i<counter; i++)
  {
    storage[i] /= (windowSize/2);
  }
  printf("Magnitude\n");

  for (i=0; i< counter; i++)
  {
    magnitude[i] = cabs(storage[i]);//(conj(storage[i])*(storage[i]));//*scale; //cabs(storage[counter]);
    //frequencies[i] =  i*(sample_freq/windowSize);
  }

  /*
  FILE *ofile= fopen("fourier.dat", "w");

  printf("Total number of elements computed %d\n", counter);
  for(i=0; i<counter ; i++)
  {
     fprintf(ofile,"%.8f\n", magnitude[i]);

  }*/
  //exit(0);

  /*Correlation among  windows
  //multiply all the elments of one window times the correspictive elments of the other one
  //sum them all
  //divide by number of elmeents (1024?)
  */

  /*
  for ( i=0; i<n_samples ; i++)
  {
    for (j=0; j<windowSize;j++)
    {
      //here I am computing the cross-correlation for windows 100 with all the other windows
      product = storage[3*windowSize]*conj(storage[i*windowSize+j]);
      summa += product;
      //correlation_result[i*windowSize+j] = product;
    }
    correlation_result[i] =  summa/windowSize;
    //printf("Sum for the window %d is:%.2f\n", i, summa);
    //normalize
    summa = 0.0;
  }
  fftw_execute(plan_backward);

  for (i=0; i<windowSize;i++)
  {
    printf("Test %.2f\n", ifft_data[i]);
  }
  */
  //to compute the correlation do :
  //e.g. select element X of the storage
  //multiply this element by all the element in a window sample
  //sum all these values
  // divide each moltiplication by this value

  /*correlation  for (i=0; i<windowSize; i++)
  {
    //conjugate = (complex double)fft_result[i+1];
    //printf("Elements R %.2f, I %.2f\n",creal(fft_result[i]), cimag(fft_result[i]));
    //correlation = fft_result[i]*conj(fft_result[i+1])*scale;
    correlation_result[counter] = correlation ;
    counter+=1;
    //printf("%2.2f, %2.2f\n", fft_result[i][0], fft_result[i][1]);
    //printf("%2.2f\n", fft_result[i][0]);
    //save fft_result in a data structure
    //test teh correlation

    //printf("Re: %.2f  , Im: %.2f \n", fft_result[i][0], fft_result[i][1]);
    //printf("Re: %2.f  , Im: %.2f \n",fft_result[i+1][0], fft_result[i+1][1]);
  }
  */

  //sanity check

  /*printf("These are the wav_data\n");
  for (i=0; i< windowSize; i++)
  {
    printf("%.2f\n", wav_data[i]);
  }
  exit(0);*/

  /*
  for (i=0; i<windowSize; i++)
  {
    stft_data[i][0] = (double)wav_data[i];
  }
  //stft_data[0][1] = 0.0 ;

  //sanity check
  printf("This is the wav_data values\n");
  printf("%.2f\n", wav_data[0]);
  printf("This is the translation into a complex object\n");
  printf("%.2f\n", stft_data[0]);
  //for ( i=0; i<windowSize;i++)
  //{
  //  printf("%.2f\n", stft_data[i]);
  //}
  //perform the fourier transform4

  printf("Fourier transform call...\n");
  fftw_execute(plan_forward);
  printf("Fourier transform done\n");
  printf("%.2f\n", fft_result[0]);
  //for (i=0; i<windowSize/2; i++)
  //{
  //  printf("%.2f\n", fft_result[i]);
  //}
  */
  //clean up the memory


  fftw_destroy_plan(plan_forward);
  fftw_free(stft_data);
  fftw_free(fft_result);
  //free(hamming_result);

  return magnitude;

}
