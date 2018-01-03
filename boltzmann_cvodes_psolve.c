#include "boltzmann_structs.h"
#include "boltzmann_cvodes_bsolve.h"
#include "boltzmann_cvodes_fsolve.h"
#include "boltzmann_cvodes_headers.h"
int boltzmann_cvodes_psolve(double t,
			    N_Vector y,
			    N_Vector fy,
			    N_Vector r,
			    N_Vector z,
			    double gamma,
			    double delta,
			    int lr,
			    void *user_data,
			    N_Vector tmp) {
  /*
    Cvodes preconditioner routine:
    Solve Pz = r
    Called by: CVode, boltzmann_cvodes_init
    Calls:     boltzmann_cvodes_bsolve, boltzmann_cvodes_fsolve
  */
  struct state_struct *state;
  double *r_data;
  double *z_data;
  int ret_code;
  int success;
  int prec_choice;
  int padi;
  state       = (struct state_struct *)user_data;
  prec_choice = state->cvodes_prec_choice;
  r_data  = NV_DATA_S(r);
  z_data = NV_DATA_S(z);
  switch (prec_choice) {
  case 0:
    /*
      This would be the diagonal preconditioner.
    */
    break;
  case 1:
    /*
      This could be the traditional ilu solver, but it
      would have an L and a U similar to case 2, so let
      it drop through.
    */
  case 2:
    if (lr == 2) {
      /*
	P = LU, backward solve z  = U^(-1)r
      */
      success = boltzmann_cvodes_bsolve(state,r_data,z_data);
    } else if (lr == 1) {
      /*
	Forward solve z = L^(-1)r
      */
      success = boltzmann_cvodes_fsolve(state,r_data,z_data);
    }
    break;
  }
  if (success) {
    ret_code = 0;
  } else {
    ret_code  = -1;
  }
  return(ret_code);
}
