#include "system_includes.h"
#include "vec_set_constant.h"
#include "ode23tb_update_wt.h"
void ode23tb_update_wt(int normcontrol, int ny, double normy, 
		     double *y, double *wt) {
  /*
    Set the wt vector to be max of the norm of y and the threshold when
    normconrol is nonzero and
    the maximum of the the absolute value of y and threshold when 
    normcontrol is 0.
    Called by: ode23tb
    Calls:     fabs
  */
  double wta;
  int    i;
  int    padi;
  if (normcontrol) {
    if (normy > wt[0]) {
      vec_set_constant(ny,wt,normy);
    }
  } else {
    for (i=0;i<ny;i++) {
      wta = fabs(y[i]);
      if (wta > wt[i])  {
	wt[i] = wta;
      }
    }
  }
}
