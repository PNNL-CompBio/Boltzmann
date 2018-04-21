#include "system_includes.h"
#include "ode23tb_normyp_o_wt.h"
#include "ode23tb_limit_h.h"
#include "ode23tb_init_h.h"
void ode23tb_init_h(int normcontrol, 
		    int ny, 
		    double normyp,
		    double recip_70p,
		    double recip_cube_root_rtol, 
		    double htspan, 
		    double hmin, 
		    double hmax, 
		    double tdir,
		    double *yp,
		    double *wt,
		    double *h_p,
		    double *absh_p) {
/*
  Compute an initial step from the  yp, vector, its 2-normb, the 
  htspan, hmin, hmax, tdir, recip70p, and recip_cube_root_rtol scalars.
  Called by: ode23tb
*/
  double yp_o_wt;
  double normyp_o_wt;
  double rh;
  double h;
  double absh;

  normyp_o_wt = ode23tb_normyp_o_wt(norm_control,ny,normyp,yp,wt);
  rh = recip_70p * normyp_o_wt * recip_cube_root_rtol;
  /*
    absh = min(hmax, htspan);
  */
  *absh_p = htspan;
  ode23tb_limit_h(hmax,hmin,rh,tdir,absh_p,h_p);
}
