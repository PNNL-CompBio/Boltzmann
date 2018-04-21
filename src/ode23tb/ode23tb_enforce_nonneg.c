#include "system_includes.h"
#include "blas.h"
#include "ode23tb_enforce_nonneg.h"
void ode23tb_enforce_nonneg(int ny, int normcontrol, double *ynew,
			    double *znew, double *normynew_p) {
  /*
    Set negative values of ynew and corresponding znew values to 0.
    If any negative values were found recopute the norm of ynew.
    Called by: ode23tb
    Calls:     dnrm2, idamax, fabs
  */
  double normynew;
  int i;
  int nnreset_znew;
  int inc1;
  int padi;
  nnreset_znew = 0;
  inc1 = 1;
  for (i=0;i<ny;i++) {
    if (ynew[i] < 0.0) {
      ynew[i] = 0.0;
      znew[i] = 0.0;
      nnreset_znew = 1;
    }
  }
  if (nnreset_znew) {
    if (normcontrol) {
      normynew = dnrm2_(&ny,ynew,&inc1);
    } else {
      normynew = fabs(ynew[idamax_(&ny,ynew,&inc1)-1]);
    }
    *normynew_p = normynew;
  }
}
