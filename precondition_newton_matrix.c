#include "boltzmann_structs.h"
#include "iluvf.h"
#include "precondition_newton_matrix.h"
int precondition_newton_matrix(struct state_struct *state) {
  /*
    Compute an approximate factorization of the newton matrix
    assumed to be stored in the miter_m, miter_im, miter_jm fields
    of the cvodes_params field of state.
    choice of preconditioner is controlled by cvodes_prec_choice field of
    state.
    Produces lower and upper triangluar factors in prec_l, prec_il, pred_jl,
    prec_u, prec_iu, and prec_ju arrays (also fields of cvodes_params).
    Called by: boltzmann_cvodes_psetup
  */
  int choice;
  int success;

  choice        = state->cvodes_prec_choice;
  success       = 1;  
  switch (choice) {
  case 0:
    /*
      nopreconditioning.
    */
    break;
  case 1:
    /*
      This ought to correspond to a diagonal preconditioner.
    */
    success = 0;
    break;
  case 2:
    /*
      This might be the traditional ILU preconditioner.
    */
    success = 0;
    break;
  case 3:
    /*
      This will be an ILU preconditioenr with value based fill.
      Number of elements allowed per factor per fill row beyond
      original number is the fill level, so a fill level of zero
      will have the same number of elements as the newtion iteraton
      matrix M.
    */
    success  = iluvf(state);
    break;
  }
  return(success);
}
