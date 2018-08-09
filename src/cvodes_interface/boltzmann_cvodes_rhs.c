#include "boltzmann_structs.h"
#include "boltzmann_cvodes_headers.h"
#include "gradient.h"
#include "boltzmann_cvodes_rhs.h"
int boltzmann_cvodes_rhs(double t, N_Vector y, N_Vector y_dot,
			 void *user_data) {
  /*
    Called by: Cvode, boltzmann_cvodes
    Calls:     gradient
  */
  struct state_struct *state;
  
  int choice;
  int ret_code;
  double *concs;
  double *flux;
  state = (struct state_struct *)user_data;
  choice = state->gradient_choice;
  concs = NV_DATA_S(y);
  flux  = NV_DATA_S(y_dot);
  ret_code = 0;
  gradient(state, concs, flux, choice);
  return(ret_code);
}
