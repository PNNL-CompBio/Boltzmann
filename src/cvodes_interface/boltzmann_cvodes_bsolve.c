#include "boltzmann_structs.h"
#include "boltzmann_cvodes_headers.h"
#include "cvodes_params_struct.h"
#include "boltzmann_cvodes_bsolve.h"
int boltzmann_cvodes_bsolve(struct state_struct *state, double *r, double *z) {
  /*
    Perform the backward solve of Uz = r
    Called by: boltzman_cvodes_psolve
  */
  struct cvodes_params_struct *cvodes_params;
  double *u;
  double *recip_diag_u;
  double sum;
  int    *iu;
  int    *ju;

  int    ny;
  int    i;

  int    j;
  int    k;
  
  int    success;
  int    padi;

  success       = 1;
  ny            = state->nunique_molecules;
  cvodes_params = state->cvodes_params;
  u             = cvodes_params->prec_u;
  iu            = cvodes_params->prec_iu;
  ju            = cvodes_params->prec_ju;
  recip_diag_u  = cvodes_params->recip_diag_u;
 
  z[ny-1] = r[ny-1]*recip_diag_u[ny-1];
  for (i=ny-2;i>=0;i--) {
    sum = r[i];
    for (j=iu[i];j<iu[i+1];j++) {
      k = ju[j];
      sum -= z[k] * u[j];
    }
    z[i] = sum * recip_diag_u[i];
  }
  return(success);
}
