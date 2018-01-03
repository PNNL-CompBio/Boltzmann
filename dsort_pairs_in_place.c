#include "system_includes.h"
#include "dsort_pairs_in_place.h"
void dsort_pairs_in_place(int n,double *list) {
  double temp;
  int i;
  int nm1;
  nm1 = n - 1;
  for (i=0;i<nm1;i+=2) {
    if (list[i] > list[i+1]) {
      temp = list[i];
      list[i] = list[i+1];
      list[i+1] = temp;
    }
  }
}

