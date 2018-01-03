#include "boltzmann_structs.h"
#include "vec_mul.h"
void vec_mul(int *n_p,double *a, double *b, double *c) {
  /*
    a <- b .* c
  */
  int i;
  int n;
  n = *n_p;
  for (i=0;i<n;i++) {
    a[i] = b[i]*c[i];
  }
}
