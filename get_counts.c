#include "system_includes.h"
#include "get_counts.h"
void get_counts(int n, double *concs, double *conc_to_count, double *counts) {
  /*
    Convert concentrations to counts.
    Called by ode23tb to enable likelihood computations.
  */
  int i;
  int padi;
  for (i=0;i<n;i++) {
    counts[i] = concs[i] * conc_to_count[i];
  }
}
