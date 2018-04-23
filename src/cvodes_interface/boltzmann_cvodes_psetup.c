#include "boltzmann_structs.h"
#include "boltzmann_cvodes_headers.h"
#include "approximate_jacobian.h"
#include "build_newton_matrix.h"
#include "precondition_newton_matrix.h"
#include "boltzmann_cvodes_psetup.h"
int boltzmann_cvodes_psetup(double t, 
			    N_Vector y,
			    N_Vector fy,
			    int jok,
			    int *jcurptr,
			    double gamma,
			    void *user_data,
			    N_Vector tmp1,
			    N_Vector tmp2,
			    N_Vector tmp3){
  /*
    Called by: CVode, boltzmann_cvodes_init
    Calls:     
  */
  struct state_struct *state;
  int success;
  int retcode;
  int choice;
  double *y_data;
  double *fy_data;
  state = (struct state_struct *) user_data;
  choice = state->ode_jacobian_choice;
  y_data = NV_DATA_S(y);
  fy_data = NV_DATA_S(fy);
  success = 1;
  retcode = 0;
  *jcurptr = 0;
  if (jok == 0) {
    /*
      Recompute jacaobian and approximate newton matrix (I - gamma J) 
      preconditioner information.
      This will set state fields dfdy_a, dfdy_ia, dfdy_ja containing the
      the jacobian matrix in compressed row storage format.
      We might want to allow for returning a dense dfdy for use in 
      ode23tb - might use a ode_jacobian_style = 0 for dense, 1 for sparse.
    */
    success = approximate_jacobian(state, y_data, fy_data, t, choice);
    if (success) {
      /*
        Build preconditioner for Newton iteration matrix.
	First build Newton iteration matrix. 
	Setting fields miter_a, miter_ia, miter_ja the compressed row storage
	format of M = (I - gamma*J)
      */
      success = build_newton_matrix(state,gamma,choice);
    }
    if (success) {
      /*
	Build the preconditioner for M.
	Here we want to build a prec_l, prec_il, prec_jl, 
	prec_u, prec_iu, and prec_ju structs to hold the approximate
	factorization of M as the preconditioner. These will be 
	used by the boltzmann_cvodes_psolve routine.
      */
      success = precondition_newton_matrix(state);
    }
    if (success) {
      *jcurptr = 1;
      retcode = 0;
    } else {
      retcode = -1;
    }
  }
  return(retcode);
}
