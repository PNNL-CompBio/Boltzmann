#include "boltzmann_structs.h"

#include "update_rxn_likelihoods.h"

#include "lr_approximate_delta_concs.h"

int lr_approximate_delta_concs(struct state_struct *state, double *counts,
			       double *forward_rxn_likelihoods,
			       double *reverse_rxn_likelihoods, 
			       double *flux, double multiplier,
			       int base_rxn, int choice) {
  /*
    Compute approximations to concentration changes wrt time, 
    Based on likelihood ratios and assumption that all
    reverse reaction likelihoods are the same as the base reaction.
    Get reference from Bill Cannon
    Called by: ode23tb, num_jac_col, ode_it_solve
    Calls:     update_rxn_likelihoods

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
         
    choice                      IOI   Not used by this routine.

    Note that multiplier is K_f(base_rxn_reaction)*(product of reactant 
    concentrations in base reaction).
	    molecule = (struct molecule_struct *)&sorted_molecules[si];
  */
  struct molecule_struct *molecules;
  struct molecule_struct *molecule;
  struct molecules_matrix_struct *molecules_matrix;
  struct rxn_matrix_struct *rxn_matrix;
  double frb;
  double lrb;
  double recip_frb;
  double forward;
  double backward;
  double *product_term;
  double *reactant_term;
  double  pt;
  double  rt;
  int64_t *molecules_ptrs;
  int64_t *rxn_indices;
  int64_t *coefficients;
  int64_t *rxn_ptrs;
  int64_t *molecule_indices;
  int64_t *rcoefficients;
  int num_species;
  int num_rxns;
  int rxn;

  int i;
  int j;

  int mi;
  int success;

  FILE *lfp;
  FILE *efp;
  /*
#define DBG 1
  */
  /*
    Check that base_rxn is in range.
  */
  success = 1;
  num_rxns = state->number_reactions;
  num_species = state->nunique_molecules;
  molecules   = state->sorted_molecules;
  molecules_matrix = state->molecules_matrix;
  molecules_ptrs   = molecules_matrix->molecules_ptrs;
  rxn_indices      = molecules_matrix->rxn_indices;
  coefficients     = molecules_matrix->coefficients;
  rxn_matrix       = state->reactions_matrix;
  rxn_ptrs         = rxn_matrix->rxn_ptrs;
  molecule_indices = rxn_matrix->molecules_indices;
  rcoefficients    = rxn_matrix->coefficients;
  product_term     = state->product_term;
  reactant_term    = state->reactant_term;

  lfp      = state->lfp;
  if ((base_rxn < 0)  || (base_rxn >= num_rxns)) {
    success = 0;
    if (lfp) {
      fprintf(lfp,"lr_approximate_delta_concs: Error, base_rxn = %d is not in [0:%d)\n",
	      base_rxn,num_rxns);
      fflush(lfp);
    }
  }
  if (success) {
    success = update_rxn_likelihoods(state,counts,forward_rxn_likelihoods,
				     reverse_rxn_likelihoods);
  }
  /*
    Fill the reactant_terms and product_terms vectors for each reaction.
    We really just need to know whether or not there are nonzero concentrations
    for species on each side of each reaction, so as not to use our
    flux approximations that assume presence of species.
  */
  for (i=0;i<num_rxns;i++) {
    pt = 1.0;
    rt = 1.0;
    for (j=rxn_ptrs[i];j<rxn_ptrs[i+1];j++) {
      mi = molecule_indices[j];
      if (rcoefficients[j] < 0) {
	rt = rt * counts[mi];
      } else {
	if (rcoefficients[j] > 0) {
	  pt = pt * counts[mi];
	}
      }
    }
    product_term[i] = pt;
    reactant_term[i] = rt;
  }
  if (success) {
    frb = forward_rxn_likelihoods[base_rxn];
    lrb = reverse_rxn_likelihoods[base_rxn];
    if (frb != 0.0) {
      recip_frb = 1.0/frb;
    } else {
      success = 0;
      if (lfp) {
	fprintf(lfp,"lr_approximate_delta_concs: Error, reaction likelihood for base_rxn reaction %d is 0.\n",base_rxn);
	fflush(lfp);
      }
    }
  }
  if (success) {
    molecule = molecules;
    for (i=0;i<num_species;i++) {
      forward = 0.0;
      backward = 0.0;
      if (molecule->variable == 1) {
	for (j=molecules_ptrs[i];j<molecules_ptrs[i+1];j++) {
	  rxn = rxn_indices[j];
	  if (coefficients[j] < 0) {
	    if (reactant_term[rxn] > 0.0) {
	      forward -= forward_rxn_likelihoods[rxn];
	    }
	    if (product_term[rxn] > 0.0) {
	      backward += 1.0;
	    }
	  } else {
	    if (coefficients[j] > 0) {
	      if (reactant_term[rxn] > 0.0) {
		forward += forward_rxn_likelihoods[rxn];
	      }
	      if (product_term[rxn] > 0.0) {
		backward -= 1.0;
	      }
	    }
	  }
	} /* end for(j...) */
	flux[i] = multiplier * ((forward * recip_frb) + (backward * lrb));
      } else {
	flux[i] = 0.0;
      }
      molecule += 1; /* Caution address arithmetic heare. */
    } /* end for (i...) */
  } /* end if success */
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"Mol_index\t   count   \t    flux\n");
    for (i=0;i<num_species;i++) {
      fprintf(lfp,"%d\t%le\t%le\n",
	      i,counts[i],flux[i]);
    }
    fprintf(lfp,"rxn_n\tfrwrd_lklhd\trvrs_lklhd\n");
    for (i=0;i<num_rxns;i++) {
      fprintf(lfp,"%d\t%le\t%le\n",i,forward_rxn_likelihoods[i],
	      reverse_rxn_likelihoods[i]);
    }
    fflush(lfp);
  }
#endif
  return (success);
}
