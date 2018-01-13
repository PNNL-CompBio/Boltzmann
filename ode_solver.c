#include "boltzmann_structs.h"
#include "ode23tb.h"
#include "boltzmann_cvodes.h"
#include "compute_dfdke_dfdmu0.h"
#include "ode_solver.h"
int ode_solver (struct state_struct *state, double *concs, int choice) {
  /*
    Called by: deq_run
    Calls:     ode23tb
               boltzmann_cvodes,
               
               concs is the vector of molecule concentrations, 
	       the concs_to_counts field in state can be used
	       to convert those concentrations to counts. 
	       We may change that to be a concentrations vector and package
	       all of the other arguments into a piece of the state_struct
	       for future convenience.

               htry, nonnegative, normcontrol, and print_concs are
	       all input variables to ode23tb the default ode solver.

	       choice is meant to be the indicator of which solver
	       to use, right now we only offer ode23tb.

	       Other solvers could pass their arguments, if different
	       through state. 


	       return value is 1 for success, 0 for failure.
  */
  int success;
  int local_choice;
  int print_output;
  FILE *lfp;
  FILE *efp;
  lfp     = state->lfp;
  print_output = (int)state->print_output;
  success = 1;
  local_choice = choice;
  switch (local_choice) {
  case 0:
    success = ode23tb(state,concs);
    break;
  case 1:
    success = boltzmann_cvodes(state,concs);
    break;
  default:
    if (lfp) {
      fprintf(lfp,"ode_solver: invalid ode_solver_choice, using default\n");
      fflush(lfp);
      state->ode_solver_choice = 0;
      success = ode23tb(state,concs);
    }
  } /* end switch(local_choice) */
  if (print_output) {
    success = compute_dfdke_dfdmu0(state,concs);
  }
  return(success);
}
