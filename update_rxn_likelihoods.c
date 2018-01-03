#include "boltzmann_structs.h"
#include "rxn_likelihoods.h"
#include "update_rxn_likelihoods.h"
int update_rxn_likelihoods(struct state_struct *state, double *counts,
			   double *forward_rxn_likelihood,
			   double *reverse_rxn_likelihood) {
  /*
    Update the forward_rxn_likelihood, and reverse_rxn_likelihood
    vectors input counts field.

    Called by: deq_run, lr_approximate_delta_concs
    Calls      rxn_likelihoods

    Arguments:
    Name        TMF       Description
    state       G*I       pointer to the boltzmann state structure.
                          Used fields are :
                             number_reactions
			     reactions and its subfields
    counts      D*I       Vector of length nunique_molecules with 
                          molecule counts for each species.
    forward_rxn_likelihood
                D*O       Vector of length number_reactions of the
		          forward reaction likelihoods

    reverse_rxn_likelihood
                D*O       Vector of length number_reactions of the
		          reverse reaction likelihoods
  */
  int success;
  int forward;

  int reverse;
  int padi;
  
  success       		   = 1;
  forward       		   = 1;
  reverse       		   = -1;
  /*
  */
  rxn_likelihoods(counts,forward_rxn_likelihood,
		  state,forward);
  rxn_likelihoods(counts,reverse_rxn_likelihood,
		  state,reverse);
  return(success);
}
