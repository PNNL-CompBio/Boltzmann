#include "system_includes.h"
#include "vec_set_constant.h"
void vec_set_constant(int ny,double *y, double value) {
  /*
    Set all elements of a vector y of length ny to a specified value.
    Called by: ode23tb, ode23tb_init_wt, ode23tb_update_wt
  */
  int i;
  int padi;
  for (i=0;i<ny;i++) {
    y[i] = value;
  }
}
