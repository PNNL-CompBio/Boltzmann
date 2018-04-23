#include "boltzmann_structs.h"
#include "boltzmann_cvodes_headers.h"
#include "cvodes_params_struct.h"
#include "lr8_approximate_ys0.h"
#include "approximate_ys0.h"
int approximate_ys0(struct state_struct *state, double *concs) {
  /*
    Fill ys0v a nunique_molecules * number_reactins (ny x ns) matrix
    with initial values for the sensitivity problem.
    ys0v]i,j]] dyi/dkej  i = 0:ny-1, j=0:nr-1
    Called by: boltzmann_cvodes_init
    Calls:     lr8_approximate_ys0
  */
  struct cvodes_params_struct *cvodes_params;
  double *ys0v;
  int    success;
  int    delta_concs_choice;
  FILE   *lfp;
  FILE   *efp;
  success = 1;
  cvodes_params      = state->cvodes_params;
  delta_concs_choice = state->delta_concs_choice;
  lfp                = state->lfp;
  ys0v               = cvodes_params->ys0v;
  switch (delta_concs_choice) {
  default:
    success = 0;
    if (lfp) {
      fprintf(lfp,"approximate_ys0: No approximation for delta_concs_choice = %d\n",delta_concs_choice);
    }
    break;
  case 8:
    success = lr8_approximate_ys0(state,ys0v,concs);
    break;
  }
  return(success);
}
