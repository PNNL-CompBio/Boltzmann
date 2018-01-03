#include "boltzmann_structs.h"

#include "update_rxn_likelihoods.h"
#include "lr_approximate_delta_concs.h"
#include "lr1_approximate_delta_concs.h"
#include "lr2_approximate_delta_concs.h"
#include "lr3_approximate_delta_concs.h"
#include "ce_approximate_delta_concs.h"

#include "approximate_delta_concs.h"

int approximate_delta_concs(struct state_struct *state, double *counts,
			    double *forward_rxn_likelihoods,
			    double *reverse_rxn_likelihoods, 
			    double *flux, double flux_scaling,
			    int base_rxn, int choice) {
  /*
    Compute approximations to concentartion changes wrt time
    0 for lr_approximate_delta_concs, based on likelihood rations.
    1 for lr1_approximate_delta_concs, based on likelihood rations.
    2 for lr2_approximate_delta_concs, based on likelihood rations.
    3 for lr3_approximate_delta_concs, based on likelihood rations.

    42 for ce_approximate_delta_concs, for debugging only with
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

    flux_scaling                  D0I   forward rate constant for base
                                      reaction multplied by base reaction
				      reactant concentration prodeuct.
				      
    base_rxn                    I0I   Base reaction number (usually 0)

    choice                      I0I   0 for lr_approximate_delta_concs
                                      2 for ce_approximate_delta_concs

    Note that flux_scaling is K_f(base_rxn_reaction)*(product of reactant 
    concentrations in base reaction).
	    molecule = (struct molecule_struct *)&sorted_molecules[si];
  */
  double alt_flux_scaling;
  int success;
  int padi;
  FILE *lfp;
  FILE *efp;
  lfp = state->lfp;
  success = 1;
  switch(choice) {
  case 0:
    /*
      Default Likelihood ratio approximation.
    */
    success = lr_approximate_delta_concs(state, 
					 counts,
					 forward_rxn_likelihoods,
					 reverse_rxn_likelihoods, 
					 flux, 
					 flux_scaling,
					 base_rxn, 
					 choice);
    break;
  case 1:
    success = lr1_approximate_delta_concs(state, 
					 counts,
					 forward_rxn_likelihoods,
					 reverse_rxn_likelihoods, 
					 flux, 
					 flux_scaling,
					 base_rxn, 
					 choice);
    break;
  case 2:
    success = lr2_approximate_delta_concs(state, 
					 counts,
					 forward_rxn_likelihoods,
					 reverse_rxn_likelihoods, 
					 flux, 
					 flux_scaling,
					 base_rxn, 
					 choice);
    break;
  case 3:
    success = lr3_approximate_delta_concs(state, 
					 counts,
					 forward_rxn_likelihoods,
					 reverse_rxn_likelihoods, 
					 flux, 
					 flux_scaling,
					 base_rxn, 
					 choice);
    break;
  case 42:
    /*
      Use debugging flavor with kinetic rate constants for coupledenzyme.in
    */
    if (state->flux_scaling == 0.0) {
      alt_flux_scaling = 1.0;
    } else {
      alt_flux_scaling = flux_scaling;
    }
    success = ce_approximate_delta_concs(state, 
					 counts,
					 forward_rxn_likelihoods,
					 reverse_rxn_likelihoods, 
					 flux, 
					 alt_flux_scaling,
					 base_rxn, 
					 choice);
    break;
  default :
    success = 0;
    if (lfp) {
      fprintf(lfp,"approximate_delta_concs: invalid delta_concs_choice was %d\n",choice);
      fflush(lfp);
    }
  } /* end switch (choice) */
  return(success);
}
