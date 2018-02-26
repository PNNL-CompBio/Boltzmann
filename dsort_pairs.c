#include "system_includes.h"
#include "dsort_pairs.h"
void dsort_pairs(int n,double *list, double *sorted_pairs) {
  /*
    Sort pairs of doubles into ascending order, with result
    sorted pairs put into sorted_pairs.
  */
  int i;
  int nm1;
  int n_odd;
  int padi;
  nm1 = n - 1;
  for (i=0;i<nm1;i+=2) {
    if (list[i] > list[i+1]) {
      sorted_pairs[i] = list[i+1];
      sorted_pairs[i+1] = list[i];
    } else {
      sorted_pairs[i] = list[i];
      sorted_pairs[i+1] = list[i+1];
    }
  }
  n_odd = n & 1;
  if (n_odd) {
    sorted_pairs[n-1] = list[n-1];
  }
}
