#include "system_includes.h"
#include "ode23tb_max_abs_ratio.h"
double ode23tb_max_abs_ratio(int n,double *a, double *b) {
  /*
    Return the maximum of (abs(a[i]/bi]) i in [0:n-1].
    Called by: ode23tb
    Calls:     fabs
  */
  double largest_r;
  double r;
  int i;
  int padi;
  largest_r = 0.0;
  for (i=0;i<n;i++) {
    r = fabs(a[i]/b[i]);
    if (r > largest_r) {
      largest_r = r;
    }
  }
  return(largest_r);
}
