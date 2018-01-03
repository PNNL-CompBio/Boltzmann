#include "boltzmann_structs.h"

#include "compute_flux_scaling.h"
#include "update_rxn_likelihoods.h"

#include "lr5_approximate_delta_concs.h"

int lr5_approximate_delta_concs(struct state_struct *state,
				double *concs,
				double *flux,
				int choice) {
  /*
    Compute approximations to concentration changes wrt time, 
    Based on likelihood ratios and assumption that all
    reverse reaction likelihoods are the same as the base reaction.
    Get reference from Bill Cannon
    Called by: approximate_delta_concs
    Calls:     update_rxn_likelihoods

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
  double  *forward_rxn_likelihoods;
  double  *reverse_rxn_likelihoods; 
  double  *counts;
  double  *conc_to_count;
  double  flux_scaling;
  double  lr_total;
  double  recip_lr_total;
  int64_t *molecules_ptrs;
  int64_t *rxn_indices;
  int64_t *coefficients;
  int64_t *rxn_ptrs;
  int64_t *molecule_indices;
  int64_t *rcoefficients;
  int num_species;
  int num_rxns;
  int rxn;
  int base_rxn;

  int i;
  int j;

  int success;
  int padi;

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
  base_rxn         = (int)state->base_reaction;
  counts           = state->ode_counts;
  conc_to_count    = state->conc_to_count;
  forward_rxn_likelihoods = state->ode_forward_lklhds;
  reverse_rxn_likelihoods = state->ode_reverse_lklhds;
  flux_scaling     = compute_flux_scaling(state,concs);

  lfp      = state->lfp;
  for (i=0;i<num_species;i++) {
    counts[i] = concs[i] * conc_to_count[i];
  }
  success = update_rxn_likelihoods(state,counts,forward_rxn_likelihoods,
				   reverse_rxn_likelihoods);
  /*    We really just need to know whether or not there are nonzero concentrations
    for species on each side of each reaction, so as not to use our
    flux approximations that assume presence of species.
  */
  lr_total = 0.0;
  for (j=0;j<num_rxns;j++) {
    lr_total += forward_rxn_likelihoods[j];
    lr_total += reverse_rxn_likelihoods[j];
  }
  if (lr_total != 0) {
    recip_lr_total = 1.0/lr_total;
  } else {
    recip_lr_total = 0.0;
    if (lfp) {
      fprintf(lfp,"lr5_approximate_delta_concs likelihoods sum is zero.\n");
      fflush(lfp);
    }
    success = 0;
  }
  if (success) {
    molecule = molecules;
    for (i=0;i<num_species;i++) {
      flux[i] = 0.0;
      if (molecule->variable == 1) {
	for (j=molecules_ptrs[i];j<molecules_ptrs[i+1];j++) {
	  rxn = rxn_indices[j];
	  flux[i] += (forward_rxn_likelihoods[rxn] -
		      reverse_rxn_likelihoods[rxn]) * coefficients[j];
	} /* end for(j...) */
	flux[i] = flux[i] * recip_lr_total * flux_scaling;
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
