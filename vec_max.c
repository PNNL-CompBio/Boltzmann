#include "boltzmann_structs.h"
#include "vec_max.h"
void vec_max(int *n_p,double *a, double *b, double *c) {
  /*
    a <- max(b,c)
  */
  int i;
  int n;
  n = *n_p;
  for (i=0;i<n;i++) {
    a[i] = (b[i] > c[i]) ? b[i] : c[i];
  }
}
