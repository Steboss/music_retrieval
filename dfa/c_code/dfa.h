#ifndef DFA_H_
#define DFA_H_

int floor_div(int a, int b);
float vector_avg(float *X, int n_elems);
float* cumsum_averaged(float* X, int n_elems, float avg, float* Xcumsum);
void dfa(float* X, int n_elems, int* scales, int scales_len, float* fluct, float*coeffs);

#endif
