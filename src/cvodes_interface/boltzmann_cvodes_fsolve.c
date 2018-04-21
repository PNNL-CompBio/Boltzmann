#include "boltzmann_structs.h"
#include "boltzmann_cvodes_headers.h"
#include "cvodes_params_struct.h"
#include "boltzmann_cvodes_fsolve.h"
int boltzmann_cvodes_fsolve(struct state_struct *state, double *r, double *z) {
  /*
    Perform the forward solve of Lz = r where L is unit lower triangular,
    and the diagonal elements are implied not stored.
    Called by: boltzman_cvodes_psolve
  */
  struct cvodes_params_struct *cvodes_params;
  double *l;
  double sum;
  int    *il;
  int    *jl;

  int    ny;
  int    i;

  int    j;
  int    k;
  
  int    success;
  int    padi;

  success       = 1;
  ny            = state->nunique_molecules;
  cvodes_params = state->cvodes_params;
  l             = cvodes_params->prec_l;
  il            = cvodes_params->prec_il;
  jl            = cvodes_params->prec_jl;
  z[0] = r[0];
  for (i=1;i<ny;i++) {
    sum = r[i];
    for (j=il[i];j<il[i+1];j++) {
      k = jl[j];
      sum -= z[k] * l[j];
    }
    z[i] = sum;
  }
  return(success);
}
