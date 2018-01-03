#include "boltzmann_structs.h"
#include "boltzmann_flatten_vgrng_state.h"
#include "boltzmann_load_agent_data.h"
int boltzmann_load_agent_data(struct state_struct *state, void *agent_data) {
  /*
    Transfer agent data into the state struct. This routine is to 
    facilitate the agent_data interface to biocellion.
    Called by: boltzmann_run
    Calls      boltzmann_flatten_vgrng_state
    
  */
  double *dagent_data;
  double *input_counts;
  double *current_counts;
  int64_t *lagent_data;
  int64_t *fvgrng_state;
  int64_t *fvgrng2_state;
  int64_t nunique_molecules;
  int64_t i;
  int success;
  int direction;
  int vgrng_state_length;
  int padi;

  success           = 1;
  nunique_molecules = state->nunique_molecules;
  lagent_data       = (int64_t*)agent_data;
  dagent_data       = (double*) agent_data;
  fvgrng_state      = lagent_data;
  fvgrng2_state     = (int64_t*)&fvgrng_state[16];
  direction = 1;
  success = boltzmann_flatten_vgrng_state(fvgrng_state,state->vgrng_state,
					  direction,&vgrng_state_length);
  if (success) {
    success = boltzmann_flatten_vgrng_state(fvgrng2_state,state->vgrng2_state,
					    direction,&vgrng_state_length);
  }
  if (success) {
    /*
      Load the input counts from the agent data.
    */
    input_counts = (double*)&dagent_data[32];
    current_counts = state->current_counts;
    for (i=0;i<nunique_molecules;i++) {
      current_counts[i] = input_counts[i];
    }
  }
  return(success);
}
