#include "boltzmann_structs.h"
#include "boltzmann_save_agent_data.h"
#include "boltzmann_build_agent_data_block.h"
int boltzmann_build_agent_data_block(struct state_struct *state, 
				     void **agent_data_p) {
  /*
    Allocate and fill an agent_data block from a boltzmann state_struct.
    Also sets the agent_data_length field of state.

    Called by: boltzmann_save_state, boltzmann_load_state
    Calls:     boltzmann_save_agent_data
  */
  void *agent_data;
  int64_t agent_data_length;
  int64_t nunique_molecules;
  int64_t number_reactions;
  int64_t one_l;
  int  success;
  int  padi;
  FILE *lfp;
  FILE *efp;
  success           = 1;
  agent_data        = NULL;
  nunique_molecules = state->nunique_molecules;
  number_reactions  = state->number_reactions;
  lfp               = state->lfp;
  one_l             = (int64_t)1;
  agent_data_length = (32 + (3*nunique_molecules) + number_reactions) * ((int64_t)8);
  agent_data = (void*)calloc(one_l,agent_data_length);
  if (agent_data == NULL) {
    success = 0;
    if (lfp) {
      fprintf(lfp,"boltzmann_build_agent_data_block: Error could not allocate %ld bytes for agent_state\n",agent_data_length);
      fflush(lfp);
    }
  } else {
    state->agent_data_length = agent_data_length;
    success = boltzmann_save_agent_data(state,agent_data);
    if (success) {
      *agent_data_p = agent_data;
    } else {
      free(agent_data);
      *agent_data_p = NULL;
    }
  }
  return(success);
}
