#include "boltzmann_structs.h"

#include "update_regulation.h"
#include "update_regulations.h"
int update_regulations(struct state_struct * state) {
  /*
    Update all the regulated activities.
    Called by: candidate_rxn, 
               compute_delta_g_forward_entropy_free_energy
    Calls:     update_regulation
  */
  int number_reactions;
  int rxn;
  int success;
  int padi;
  number_reactions = (int)state->number_reactions;
  success          = 1;
  for (rxn = 0;((rxn < number_reactions) && success);rxn++) {
    success= update_regulation(state,rxn);
  }
  return(success);
}
