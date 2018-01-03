#include "system_includes.h"
#include "pairwise_sum.h"
double pairwise_sum(int n, double *v, double *work) {
  /*
    Pairwise sum elements in vector v in order to reduce rounding error.
    n is length of vector v to be summed. Work is vector of length n,
    to keep partials. Called by stable_add.
  */
  double *a;
  double *b;
  double *c;
  double sum;
  int i;
  int j;
  int m;
  int cpos;
  m = n;
  a = v;
  b = work;
  cpos = (m+1) >> 1;
  c = &b[cpos];
  while (m > 1) {
    j=0;
    for (i=0;i<m-1;i+=2) {
      b[j] = a[i] + a[i+1];
      j++;
    }
    if (m & 1) {
      b[j] = a[m-1];
      j++;
    }
    m = j;
    a = b;
    b = c;
    c = a;
  }
  sum = a[0];
  return(sum);
}
