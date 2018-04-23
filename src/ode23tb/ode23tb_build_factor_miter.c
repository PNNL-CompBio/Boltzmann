#include "system_includes.h"
#include "blas.h"
#include "lapack.h"
#include "ode23tb_build_factor_miter.h"
int  ode23tb_build_factor_miter(int ny,
				int nysq,
				double d,
				double h,
				double *dfdy,
				double *miter,
				int    *ipivot,
				int    *info_p,
				FILE   *lfp) {
  /*
    Build the iteration matrix miter
    This version assumes the mass conservation matrix is the identity.
    Called by: ode23tb
    Calls:     dscal_,daxpy_,dgetrf_,
               fprintf, fflush
  */
  double dzero;
  double mdh;
  int info;
  int inc1;
  int success;
  int i;
  success = 1;
  inc1    = 1;
  dzero   = 0.0;
  /*
    Set miter to be all 0.
  */
  dscal_(&nysq,&dzero,miter,&inc1);
  /*
    Set miter to be the identity matrix.
  */
  for (i=0;i<nysq;i += (ny+1)) {
    miter[i] = 1.0;
  }
  /*
    set miter to be I - (d*h)*dfdy
  */
  mdh = 0.0 - (d*h);
  daxpy_(&nysq,&mdh,dfdy,&inc1,miter,&inc1);
  /*
    Now we need to factor miter with say dgetrf (lapack routines)
    overwrites miter with its factorization, generating
    a pivot vector, ipivot,
    as arguments to ode_it_solve, 
  */
  info = 0;
  dgetrf_(&ny,&ny,miter,&ny,ipivot,&info);
  if (info != 0) {
    if (lfp) {
      fprintf(lfp,"ode23tb_build_factor_miter: "
	      "dgetrf_ call failed with info = %d\n",
	      info);
      fflush(lfp);
    }
    success = 0;
  }
  return(success);
}
