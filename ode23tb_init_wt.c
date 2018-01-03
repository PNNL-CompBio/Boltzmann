#include "system_includes.h"
#include "vec_set_constant.h"
#include "ode23tb_init_wt.h"
void ode23tb_init_wt(int normcontrol, int ny, double normy, double threshold,
		     double *y, double *wt) {
  /*
    Set the wt vector to be max of the norm of y and the threshold when
    normconrol is nonzero and
    the maximum of the the absolute value of y and threshold when 
    normcontrol is 0.
    Called by: ode23tb
    Calls:     vec_set_constant.fabs
  */
  double wt2;
  double wta;
  int    i;
  int    padi;
  if (normcontrol) {
    wt2 = normy;
    if (threshold > wt2) {
      wt2 = threshold;
    }
    vec_set_constant(ny,wt,wt2);
  } else {
    for (i=0;i<ny;i++) {
      wta = fabs(y[i]);
      if (threshold > wta) {
	wta = threshold;
      }
      wt[i] = wta;
    }
  }
}
