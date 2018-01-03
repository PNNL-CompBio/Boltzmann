#include "boltzmann_structs.h"
#include "blas.h"
#include "build_newton_matrix.h"
int  build_newton_matrix(struct state_struct *state,
			 double gamma,
			 int choice) {
  /*
    Build a sparse Newton iteration matrix from the
    sparse jacobian stored in dfdy_a, dfdy_ia, dfdy_ja, and gamma
    M = I - gamma J
    Called by: boltzmann_cvodes_psetup
  */
  struct cvodes_params_struct *cvodes_params;
  double *miter_m;
  double *dfdy_a;
  double mgamma;
  double one;
  int    *miter_im;
  int    *miter_jm;
  int    *dfdy_ia;
  int    *dfdy_ja;

  int    nnz;
  int    ny;

  int    inc1;
  int    success;

  int    a_len;
  int    ja_len;

  int    ia_len;
  int    i;

  int    j;
  int    k;

  success       = 1;
  ny            = state->nunique_molecules;
  cvodes_params = state->cvodes_params;
  dfdy_a        = cvodes_params->dfdy_a;
  dfdy_ia       = cvodes_params->dfdy_ia;
  dfdy_ja       = cvodes_params->dfdy_ja;
  miter_m       = cvodes_params->miter_m;
  miter_im      = cvodes_params->miter_im;
  miter_jm      = cvodes_params->miter_jm;
  nnz           = dfdy_ia[ny];
  one           = 1.0;
  mgamma        = 0.0 - gamma;
  a_len         = nnz << 3;
  ja_len        = nnz << 2;
  ia_len        = (ny+1) << 2;
  memcpy(miter_m,dfdy_a,a_len);
  memcpy(miter_im,dfdy_ia,ia_len);
  memcpy(miter_jm,dfdy_ja,ja_len);
  /*
    Scale copy of a in m by -gamma;.
  */
  inc1 = 1;
  dscal_(&nnz,&mgamma,miter_m,&inc1);
  /*
    Add 1 to the diagonal.
  */
  for (i=0;i<ny;i++) {
    for (k=miter_im[i];i<miter_im[i+1];k++) {
      j = miter_jm[k];
      if (j == i) {
	miter_m[k] += 1.0;
      }
    }
  }
  return(success);
}
