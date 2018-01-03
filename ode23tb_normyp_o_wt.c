#include "system_includes.h"
#include "ode23tb_normyp_o_wt.h"
double  ode23tb_normyp_o_wt(int normcontrol, int ny, double normyp,
			    double *yp, double *wt) {

  /*
    compute the 2 (normcontrol !=0) or infinity (normcontrol == 0)
    norm of yp./wt
    Called by: ode23tb, ode23tb_init_h
    Calls:     fabs
  */

  double normyp_o_wt;
  double yp_o_wt;
  int i;
  int padi;
  if (normcontrol) {
    normyp_o_wt = normyp/wt[0];
  } else {
    normyp_o_wt = 0.0;
    for (i=0;i<ny;i++) {
      yp_o_wt = fabs(yp[i])/wt[i];
      if (yp_o_wt > normyp_o_wt) {
	normyp_o_wt = yp_o_wt;
      }
    }
  }
  return(normyp_o_wt);
}
