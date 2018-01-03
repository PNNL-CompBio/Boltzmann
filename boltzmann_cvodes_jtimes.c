#include "boltzmann_structs.h"
#include "boltzmann_cvodes_headers.h"
#include "approximate_jacobian.h"
#include "boltzmann_sparse_mvp.h"
#include "boltzmann_cvodes_jtimes.h"
int boltzmann_cvodes_jtimes(N_Vector v,
			    N_Vector jv,
			    double t,
			    N_Vector y,
			    N_Vector fy,
			    void *user_data,
			    N_Vector tmp) {
  /*
    Multiply an N_Vector by the Jacobian approximation, first
    forming the jacobian.
    Called by: Cvode, boltzmann_cvodes_init
    Calls: approximate_jacovban,
           boltzman_sparse_mvp
  */
  struct state_struct *state;
  struct cvodes_params_struct *cvodes_params;
  double *v_data;
  double *jv_data;
  double *y_data;
  double *fy_data;
  double *dfdy_a;
  int *dfdy_ia;
  int *dfdy_ja;
  int choice;
  int success;
  int ny;
  int ret_code;
  state = (struct state_struct *)user_data;
  ny      = state->nunique_molecules;
  cvodes_params = state->cvodes_params;
  choice  = state->ode_jacobian_choice;
  ret_code = -1;
  if (choice == 0) {
    success = 0;
  } else {
    /*
      unpack N_Vectors.
    */
    v_data  = NV_DATA_S(v);
    jv_data = NV_DATA_S(jv);
    y_data  = NV_DATA_S(y);
    fy_data = NV_DATA_S(fy);
    ret_code = -1;
    /*
      First we need to build the jacobian with a call to approximate_jacobian.
    */
    success = approximate_jacobian(state, y_data, fy_data, t, choice);
    /*
      Then we need to apply the jacobian to v yeilding Jv.
    */
    if (success) {
      dfdy_a = cvodes_params->dfdy_a;
      dfdy_ia = cvodes_params->dfdy_ia;
      dfdy_ja = cvodes_params->dfdy_ja;
      success  = boltzmann_sparse_mvp(ny,dfdy_a, dfdy_ia, dfdy_ja, v_data,
				      jv_data);
    }
  }
  if (success) {
    ret_code = 0;
  }
  return(ret_code);
}
