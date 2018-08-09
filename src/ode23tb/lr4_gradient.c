#include "boltzmann_structs.h"

#include "compute_flux_scaling.h"
#include "update_rxn_likelihoods.h"
#include "conc_to_pow.h"

#include "lr4_gradient.h"

int lr4_gradient(struct state_struct *state, 
				double *concs,
				double *flux, 
				int choice) {
  /*
    Compute approximations to concentration changes wrt time, 
    Based on likelihood ratios and assumption that all
    reverse reaction likelihoods are the same as the base reaction.
    Get reference from Bill Cannon
    Called by: gradient
    Calls:     update_rxn_likelihoods,
               compute_flux_scaling,
	       conc_to_pow

                                TMF
    state                       *SI   Boltzmant state structure.
                                      uses number_reactions,
				           unique_moleules,
                                           molecules_matrix,
					   and lfp,
				      

    concs			D1I   molecule concentrations vector of length 
                                      nunique_moleucles

    flux                        D1O   vector of length  unique_molecules
                                      of concentration change per unit time.
				      Set by this routine.

    choice                      IOI   Not used by this routine.

  */
  struct  molecule_struct *molecules;
  struct  molecule_struct *molecule;
  struct  molecules_matrix_struct *molecules_matrix;
  struct  reactions_matrix_struct *rxn_matrix;
  double  *product_term;
  double  *reactant_term;
  double  *forward_rxn_likelihoods;
  double  *reverse_rxn_likelihoods; 
  double  *counts;
  double  *conc_to_count;
  double  *coefficients;
  double  *rcoefficients;
  double  *kq;
  double  *kqi;
  double  *skq;
  double  *skqi;
  double  *activities;
  int64_t *molecules_ptrs;
  int64_t *rxn_indices;
  int64_t *rxn_ptrs;
  int64_t *molecule_indices;
  double  flux_scaling;
  double  frb;
  double  lrb;
  double  recip_frb;
  double  forward;
  double  backward;
  double  pt;
  double  rt;
  double  klim;
  double  count_mi;
  double  factorial;
  int num_species;
  int num_rxns;
  int rxn;
  int base_rxn;

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
  rxn_indices      = molecules_matrix->reaction_indices;
  coefficients     = molecules_matrix->coefficients;
  rxn_matrix       = state->reactions_matrix;
  rxn_ptrs         = rxn_matrix->rxn_ptrs;
  molecule_indices = rxn_matrix->molecules_indices;
  rcoefficients    = rxn_matrix->coefficients;
  product_term     = state->product_term;
  reactant_term    = state->reactant_term;
  base_rxn         = (int)state->base_reaction;
  counts           = state->ode_counts;
  conc_to_count    = state->conc_to_count;
  forward_rxn_likelihoods = state->ode_forward_lklhds;
  reverse_rxn_likelihoods = state->ode_reverse_lklhds;
  kq               = state->ode_kq;
  kqi              = state->ode_kqi;
  skq              = state->ode_skq;
  skqi             = state->ode_skqi;
  activities       = state->activities;
  factorial        = 0.0;
  flux_scaling     = compute_flux_scaling(state,concs);
  frb              = 1.0;
  recip_frb        = 1.0;
  lfp      = state->lfp;
  if ((base_rxn < 0)  || (base_rxn >= num_rxns)) {
    success = 0;
    if (lfp) {
      fprintf(lfp,"lr4_gradient: Error, base_rxn = %d is not in [0:%d)\n",
	      base_rxn,num_rxns);
      fflush(lfp);
    }
  }
  if (success) {
    for (i=0;i<num_species;i++) {
      counts[i] = concs[i] * conc_to_count[i];
    }
    success = update_rxn_likelihoods(state,counts,forward_rxn_likelihoods,
				     reverse_rxn_likelihoods);
  }
  /*
    Fill the reactant_terms and product_terms vectors for each reaction.
    We really just need to know whether or not there are nonzero concentrations
    for species on each side of each reaction, so as not to use our
    flux approximations that assume presence of species.
  */
  if (success) {
    frb = forward_rxn_likelihoods[base_rxn];
    lrb = reverse_rxn_likelihoods[base_rxn];
    if (frb != 0.0) {
      recip_frb = 1.0/frb;
    } else {
      success = 0;
      if (lfp) {
	fprintf(lfp,"lr4_gradient: Error, reaction likelihood for base_rxn reaction %d is 0.\n",base_rxn);
	fflush(lfp);
      }
    }
  }
  if (success) {
    for (i=0;i<num_rxns;i++) {
      pt = 1.0;
      rt = 1.0;
      for (j=rxn_ptrs[i];j<rxn_ptrs[i+1];j++) {
        mi = molecule_indices[j];
        klim = rcoefficients[j];
        count_mi = counts[mi];
        if (klim < 0.0) {
	  klim = 0.0 - klim;
	  rt = rt * conc_to_pow(count_mi,klim,factorial);
	  /*
	    for (k=0;k<(-klim);k++) {
	    rt = rt * counts[mi];
	    }
	  */
        } else {
	  if (klim > 0.0) {
	    pt = pt * conc_to_pow(count_mi,klim,factorial);
	    /*
	      for (k=0;k<klim;k++) {
	      pt = pt * counts[mi];
	      }
	    */
	  }
        }
      }
      kq[i]  = forward_rxn_likelihoods[i] * recip_frb;
      kqi[i] = reverse_rxn_likelihoods[i] * recip_frb;
      skq[i] = kq[i] * activities[i];
      skqi[i] = kqi[i] * activities[i];
      product_term[i] = pt;
      reactant_term[i] = rt;
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
	  if (coefficients[j] < 0.0) {
	    if (reactant_term[rxn] > 0.0) {
	      forward -= forward_rxn_likelihoods[rxn];
	    }
	    if (product_term[rxn] > 0.0) {
	      backward += reverse_rxn_likelihoods[rxn];
	    }
	  } else {
	    if (coefficients[j] > 0.0) {
	      if (reactant_term[rxn] > 0.0) {
		forward += forward_rxn_likelihoods[rxn];
	      }
	      if (product_term[rxn] > 0.0) {
		backward -= reverse_rxn_likelihoods[rxn];
	      }
	    }
	  }
	} /* end for(j...) */
	flux[i] = flux_scaling * recip_frb * (forward + backward);
      } else {
	flux[i] = 0.0;
      }
      molecule += 1; /* Caution address arithmetic here. */
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
