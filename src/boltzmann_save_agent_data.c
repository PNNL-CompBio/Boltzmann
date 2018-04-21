#include "boltzmann_structs.h"
#include "boltzmann_flatten_vgrng_state.h"
#include "boltzmann_save_agent_data.h"
int boltzmann_save_agent_data(struct state_struct *state, void *agent_data) {
  /*
    Transfer agent data from  the state struct. This routine is to 
    facilitate the agent_data interface to biocellion.
    Called by: boltzmann_run, boltzmann_build_agent_data_block
    Calls      boltzmann_flatten_vgrng_state
    
  */
  double *dagent_data;
  double *current_counts;
  double *output_counts;
  double *bndry_flux_counts;
  double *agent_bndry_flux_counts;
  double *net_lklhd_bndry_flux;
  double *agent_net_lklhd_bndry_flux;
  double *net_likelihood;
  double *agent_net_likelihood;
  int64_t *lagent_data;
  int64_t *fvgrng_state;
  int64_t *fvgrng2_state;
  int64_t nunique_molecules;
  int64_t number_reactions;
  int64_t i;

  int success;
  int direction;
  int vgrng_state_length;
  int padi;

  success = 1;
  nunique_molecules = state->nunique_molecules;
  number_reactions  = state->number_reactions;
  lagent_data  = (int64_t*)agent_data;
  dagent_data  = (double*) agent_data;
  fvgrng_state = lagent_data;
  fvgrng2_state = (int64_t*)&fvgrng_state[16];
  direction = 0;
  success = boltzmann_flatten_vgrng_state(fvgrng_state,state->vgrng_state,
					  direction,&vgrng_state_length);
  if (success) {
    success = boltzmann_flatten_vgrng_state(fvgrng2_state,state->vgrng2_state,
					    direction,&vgrng_state_length);
  }
  if (success) {
  /*
    Save the output counts to the agent data.
  */
    output_counts = (double*)&dagent_data[32];
    current_counts = state->current_counts;
    for (i=0;i<nunique_molecules;i++) {
      output_counts[i] = current_counts[i];
    }
    bndry_flux_counts    = state->bndry_flux_counts;
    net_lklhd_bndry_flux = state->net_lklhd_bndry_flux;
    net_likelihood       = state->net_likelihood;
    agent_bndry_flux_counts = (double*)&output_counts[nunique_molecules];
    agent_net_lklhd_bndry_flux = (double*)&agent_bndry_flux_counts[nunique_molecules];
    agent_net_likelihood = (double*)&agent_net_lklhd_bndry_flux[nunique_molecules];
    for (i=0;i<nunique_molecules;i++) {
      agent_bndry_flux_counts[i] = bndry_flux_counts[i];
    }
    for (i=0;i<nunique_molecules;i++) {
      agent_net_lklhd_bndry_flux[i] = net_lklhd_bndry_flux[i];
    }
    for (i=0;i<number_reactions;i++) {
      agent_net_likelihood[i] = net_likelihood[i];
    }
  }
  return(success);
}
