double* linspace(double *x, double init, double fin, int N);
double mel_to_hz(double curr_mel);
double hz_to_mel(double frequency);
double* mel_frequencies(int n_mels, double fqmin, double fqmax, double* melfreqs);
void mel(int sample_freq,int windowSize, double** weights );
double* mfcc(double *wav_data, int samples, int windowSize, int hop_size,\
             double *magnitude, int sample_freq, int length, int n_mfcc, double *cython_basis);
