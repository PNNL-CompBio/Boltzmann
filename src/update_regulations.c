#include "boltzmann_structs.h"

#include "update_regulation.h"
#include "update_regulations.h"
int update_regulations(struct state_struct * state, double *counts_or_concs,
		       int count_or_conc) {
  /*
    Update all the regulated activities, based on the counts_or_concs
    vector of counts (if count_or_conc == 1) or concentrations if 
    (count_or_conc != 1).

    Called by: candidate_rxn, 
               compute_delta_g_forward_entropy_free_energy
	       lr8_gradient,
	       lr9_gradient,
	       lr10_gradient,
	       lr11_gradient,
	       lr12_gradient,
	       lr8_approximate_jacobian,
    Calls:     update_regulation
  */
  int number_reactions;
  int rxn;
  int success;
  int padi;
  number_reactions = (int)state->number_reactions;
  success          = 1;
  for (rxn = 0;((rxn < number_reactions) && success);rxn++) {
    success= update_regulation(state,rxn,counts_or_concs,count_or_conc);
  }
  return(success);
}
