#include "system_includes.h"
#include "ode23tb_nonneg_err.h"
void ode23tb_nonneg_err(int ny, 
			double *ynew,
			double rtol, 
			double ernn_scale,
			double *err_p,
			int *nnrejectstep_p) {
  /*
    Compute the norm of the negative ynew values and scale by ernn_scale,
    then compare against rtol, setting *nnrejectsetp_p to 1 if the 
    scaled negative value norm is > rtol.
    Called by: ode23tb
    Calls:     sqrt
  */
  double errnn;
  double err;
  int i;
  int nnrejectstep;
  err = *err_p;
  nnrejectstep = 0;
  errnn = 0.0;
  if (err <= rtol) {
    for (i=0;i<ny;i++) {
      if (ynew[i] < 0.0) {
	errnn += (ynew[i] * ynew[i]);
      }
    }
    errnn = sqrt(errnn) * ernn_scale;;
    if (errnn > rtol) {
      err = errnn;
      nnrejectstep = 1;
    }
  }
  *err_p = err;
  *nnrejectstep_p = nnrejectstep;
} /* end if nonnegative */

