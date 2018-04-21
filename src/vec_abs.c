#include "boltzmann_structs.h"
#include "vec_abs.h"
void vec_abs(int *n_p, double *a, double *b) {
/*
  a <- |b|
*/
  int i;
  int n;
  n = *n_p;
  for (i=0;i<n;i++) {
    a[i] = fabs(b[i]);
  }
}
