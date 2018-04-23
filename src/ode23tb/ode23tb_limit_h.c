#include "system_includes.h"
#include "ode23tb_limit_h.h"
void ode23tb_limit_h(double hmax, double hmin, double rh, double tdir,
		    double *absh_p,double *h_p) {
  /*
    Limit |h| to be in [hmin,min(hmax,1/rh)]
    Called by: ode23tb_init_h, ode23tb
  */
  double absh;
  double h;
  absh = *absh_p;
  if (absh > hmax) {
    absh = hmax;
  }
  if ((absh * rh) > 1.0) {
    absh = 1.0/rh;
  }
  /*
    absh = max(absh,hmin);
  */
  if (hmin > absh) {
    absh = hmin;
  }
  h = tdir * absh;
  *h_p = h;
  *absh_p = absh;
}
