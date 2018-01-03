#include "boltzmann_structs.h"
#include "vgrng_init.h"
#include "print_rxn_likelihoods_header.h"
#include "print_free_energy_header.h"
#include "flatten_state.h"
/*
#include "zero_solvent_coefficients.h"
*/
#include "free_boot_state2.h"

#include "run_init.h"
int run_init(struct state_struct *state, struct state_struct **flattened_state) {
/*
    Initialize the activities and the random number generator,
    pack the state into a new state variable, flattened_state,
    and free the fields of the boot state.
    
    Called by: boltzmann_init_core, 
    Calls:     vgrng_init,
               print_rxn_likelihoods_header,
	       print_free_endrgy_header,
               flatten_state,
	       zero_solvent_coefficients,
	       free_boot_state2

  */
  struct state_struct *stateq;
  struct vgrng_state_struct *vgrng_state;
  struct vgrng_state_struct *vgrng2_state;
  double *activities;
  int64_t vgrng_start;
  int64_t i;

  int success;
  int vgrng_start_steps;

  int print_output;
  int padi;
  
  FILE *lfp;
  FILE *efp;

  success = 1;
  print_output = state->print_output;
  activities = state->activities;
  if (state->use_activities == 0) {
    for (i=0;i<state->number_reactions;i++) {
      activities[i] = 1.0;
    }
  }
  /*
    Initialize the random number generators,
    setting the vgrng_state and vgrng2_state fields of
    the state structure.
  */
  vgrng_state = state->vgrng_state;
  vgrng_start_steps = 1001;
  vgrng_start= vgrng_init(vgrng_state,vgrng_start_steps);
  vgrng2_state = state->vgrng2_state;
  vgrng_start_steps = 1042;
  vgrng_start= vgrng_init(vgrng2_state,vgrng_start_steps);

  if (state->print_output) {
    state->rxn_view_hist_length = ((int64_t)(state->record_steps + state->rxn_view_freq -2)/state->rxn_view_freq) + (int64_t)1;
  } else {
    state->rxn_view_hist_length = 0;
  }
  if (print_output) {
    /*
      Print the header lines for the reaction likelihoods output file.
    */
    print_rxn_likelihoods_header(state);
    if (state->free_energy_format > (int64_t)0) {
      /*
	Print the header lines for the free energy output file.
      */
      print_free_energy_header(state);
    }
  }
  stateq  = NULL;
  state->workspace_base = NULL;
  success = flatten_state(state,&stateq);
  /*
  if (success) {
    success = zero_solvent_coefficients(stateq);
  }
  */
  *flattened_state = stateq;
  if (success) {
    success - free_boot_state2(state);
  }
  return(success);
}
