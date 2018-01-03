#include "boltzmann_structs.h"
#include "compute_net_likelihoods.h"
int compute_net_likelihoods(struct state_struct *state,
			    double *forward_rxn_likelihood,
			    double *reverse_rxn_likelihood,
			    double *net_likelihood) {
  /*
    Compute the net likelihood for each reaction.
    Uses the forward_rxn_likelihood and reverse_rxn_likelihood 
    fields in state to set the net_likelilhood field for each 
    reaction.
    Called by: ode23b
    Calls:
  */
  int64_t number_reactions;
  double  total_likelihood;
  double  recip_total_likelihood;
  int64_t i;
  int success;
  int padi;
  FILE *lfp;
  FILE *efp;
  success                = 1;
  number_reactions       = state->number_reactions;
  lfp                    = state->lfp;
  total_likelihood = 0;
  for (i=0;i<number_reactions;i++) {
    total_likelihood = total_likelihood + forward_rxn_likelihood[i] +
      reverse_rxn_likelihood[i];
    net_likelihood[i] = forward_rxn_likelihood[i] - reverse_rxn_likelihood[i];
  }
  if (total_likelihood <= 0) {
    success = 0;
    if (lfp) {
      fprintf(lfp,"compute_net_likelihoods: Error total likelihood is <= 0.\n");
      fflush(lfp);
    }
  } else {
    recip_total_likelihood = 1.0 / total_likelihood;
    for (i=0;i<number_reactions;i++) {
      net_likelihood[i] *= recip_total_likelihood;
    }
  }
  return(success);
}
