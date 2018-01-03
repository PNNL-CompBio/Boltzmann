#include "system_includes.h"
#include "dreverse_list.h"
void dreverse_list(int n, double *v) {
  /*
    Reverse an array of doubles.
    Called by: stable_add
  */
  double temp;
  int i;
  int m;
  int h;
  h = n>>1;
  m = n-1;
  for (i=0;i<h;i++) {
    temp = v[i];
    v[i] = v[m];
    v[m] = temp;
    m--;
  }
}
