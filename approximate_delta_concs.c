#include "boltzmann_structs.h"

#include "update_rxn_likelihoods.h"
#include "lr_approximate_delta_concs.h"
#include "ce_approximate_delta_concs.h"

#include "approximate_delta_concs.h"

int approximate_delta_concs(struct state_struct *state, double *counts,
			    double *forward_rxn_likelihoods,
			    double *reverse_rxn_likelihoods, 
			    double *flux, double multiplier,
			    int base_rxn, int choice) {
  /*
    Compute approximations to concentartion changes wrt time
    0 for lr_approximate_delta_concs, based on likelihood rations.

    2 for ce_approximate_delta_concs, for debugging only with
    the reaction rates for the coupledenzyme.in file

    Called by: ode23tb, num_jac_col, ode_it_solve
    Calls      lr_approximate_delta_concs,
               ce_approximate_delta_concs,
	       

                                TMF
    state                       *SI   Boltzmant state structure.
                                      uses number_reactions,
				           unique_moleules,
                                           molecules_matrix,
					   and lfp,
				      

    counts			D1I   molecule counts vector of length 
                                      nunique_moleucles
    forward_rxn_likelihoods     D1W   scratch vector of length number_reactions
    reverse_rxn_likelihoods     D1W   scratch vector of length number_reactions

    flux                        D1O   vector of length  unique_molecules
                                      of concentration change per unit time.
				      Set by this routine.

    multiplier                  D0I   forward rate constant for base
                                      reaction multplied by base reaction
				      reactant concentration prodeuct.
				      
    base_rxn                    I0I   Base reaction number (usually 0)

    choice                      I0I   0 for lr_approximate_delta_concs
                                      2 for ce_approximate_delta_concs

    Note that multiplier is K_f(base_rxn_reaction)*(product of reactant 
    concentrations in base reaction).
	    molecule = (struct molecule_struct *)&sorted_molecules[si];
  */
  int success;
  int padi;
  FILE *lfp;
  FILE *efp;
  lfp = state->lfp;
  success = 1;
  if (choice == 2) {
    /*
      Use debugging flavor with rate constants for coupledenzyme.in
    */
    success = ce_approximate_delta_concs(state, 
					 counts,
					 forward_rxn_likelihoods,
					 reverse_rxn_likelihoods, 
					 flux, 
					 multiplier,
					 base_rxn, 
					 choice);
  } else {
    /*
      Default Likelihood ratio approximation.
    */
    if (choice != 0) {
      if (lfp) {
	fprintf(lfp,"approximate_delta_concs: Invalid choice, using default.\n");
	fflush(lfp);
      }
      state->delta_concs_choice = (int64_t)0;
    }
    success = lr_approximate_delta_concs(state, 
					 counts,
					 forward_rxn_likelihoods,
					 reverse_rxn_likelihoods, 
					 flux, 
					 multiplier,
					 base_rxn, 
					 choice);
  }
  return(success);
}
