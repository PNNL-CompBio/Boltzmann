#include "boltzmann_structs.h"
#include "boltzmann_cvodes_headers.h"
#include "cvodes_params_struct.h"
#include "get_counts.h"
#include "conc_to_pow.h"
#include "update_regulations.h"
#include "lr14_approximate_delta_concs.h"

int lr14_approximate_delta_concs(struct state_struct *state, 
				double *concs,
				double *flux, 
				int choice) {
  /*
    Compute approximations to concentration changes wrt time, 
    Based on thermodynamics formulation for concentraion rate
    changes
    using counts to compute tr, tp, pt, rt instead of concs.

    Get reference from Bill Cannon
    Called by: approximate_delta_concs
    Calls:     get_counts,
               conc_to_pow,
               update_regulations,

                                TMF
    state                       *SI   Boltzmant state structure.
                                      uses number_reactions,
				           unique_moleules,
					   sorted_molecules,
                                           molecules_matrix,
					   ke, rke,
					   product_term as scratch.
					   and lfp,
				      

    concs			D1I   molecule concentrations vector of length 
                                      nunique_moleucles

    flux                        D1O   vector of length  unique_molecules
                                      of concentration change per unit time.
				      Set by this routine.

    choice                      IOI   Not used by this routine.

  */
  struct  cvodes_params_struct *cvodes_params;
  struct  molecule_struct *molecules;
  struct  molecule_struct *molecule;
  struct  compartment_struct *compartments;
  struct  compartment_struct *compartment;
  struct  molecules_matrix_struct *molecules_matrix;
  struct  reactions_matrix_struct *rxn_matrix;
  double  *activities;
  double  *forward_lklhd;
  double  *reverse_lklhd;
  double  *rfc;
  double  *ke;
  double  *rke;
  double  *counts;
  double  *conc_to_count;
  double  flux_scaling;
  double  pt;
  double  rt;
  double  tr;
  double  tp;
  double  volume;
  double  recip_volume;
  double  multiplier;
  double  keq_adj;
  double  rkeq_adj;
  double  coeff;
  double  avogadro;
  double  recip_avogadro;
  double  fluxi;
  double  conc_mi;
  double  factorial;
  double  klim;
  /*
  double  sum_coeff;
  */
  double  volume_avo;
  double  recip_volume_avo;
  /*
  double  count_mi;
  */
  double  *coefficients;
  double  *rcoefficients;
  int64_t *molecules_ptrs;
  int64_t *rxn_indices;
  int64_t *rxn_ptrs;
  int64_t *molecule_indices;

  double *coeff_sum;
  int num_species;
  int num_rxns;

  int rxn;
  int success;

  int i;
  int j;

  int mi;
  int ci;

  int use_regulation;
  int count_or_conc;

  int compute_sensitivities;
  int ode_solver_choice;

  FILE *lfp;
  FILE *efp;
  /*
#define DBG 1
  */
  /*
    Check that base_rxn is in range.
  */
  success          = 1;
  num_rxns         = state->number_reactions;
  num_species      = state->nunique_molecules;
  molecules        = state->sorted_molecules;
  compartments     = state->sorted_compartments;
  activities       = state->activities;
  forward_lklhd    = state->ode_forward_lklhds;
  reverse_lklhd    = state->ode_reverse_lklhds;
  coeff_sum        = state->coeff_sum;
  molecules_matrix = state->molecules_matrix;
  molecules_ptrs   = molecules_matrix->molecules_ptrs;
  rxn_indices      = molecules_matrix->reaction_indices;
  coefficients     = molecules_matrix->coefficients;
  rxn_matrix       = state->reactions_matrix;
  rxn_ptrs         = rxn_matrix->rxn_ptrs;
  molecule_indices = rxn_matrix->molecules_indices;
  rcoefficients    = rxn_matrix->coefficients;
  ke               = state->ke;
  rke              = state->rke;
  rfc              = state->rfc;
  counts           = state->ode_counts;
  conc_to_count    = state->conc_to_count;
  use_regulation   = state->use_regulation;
  ode_solver_choice = state->ode_solver_choice;
  compute_sensitivities = state->compute_sensitivities;
  factorial           = 0.0;
  if ((ode_solver_choice == 1) && compute_sensitivities) {
    cvodes_params = state->cvodes_params;
    ke = cvodes_params->p;
    rke = cvodes_params->rp;
    for (i=0;i<num_rxns;i++) {
      rke[i] = 1.0/ke[i];
    }
  }
  avogadro         = state->avogadro;
  recip_avogadro   = state->recip_avogadro;
  /*
  flux_scaling     = compute_flux_scaling(state,concs);
  */
  get_counts(num_species,concs,conc_to_count,counts);
  flux_scaling     = 1.0;
  lfp      = state->lfp;
  /*
    As per discusion with Bill Cannon, we want to update the activities
    if reguation is in play. So do that here.
  */
  if (use_regulation) {
    count_or_conc = 0;
    update_regulations(state,concs,count_or_conc);
  }
  /*
    Compute the reaction flux contributions for each reaction:

    rfc   = k_f * product of reactants^stoichiometric_coef -
            k_r * product of products^stoichiometric_coef,

	 where k_f = k_eq/ thermo_product of products.

	 and   k_r  = k_eq^(-1)/thermo_product of reactants,

	 and thermo_product = 
	 product( species_conc + |stoichiometric_coef|/volume)^|stoichiometric coef|.

	 Then if a molecule is in rxn i,

	 rxn i contributes 1/stoichiometric_coef * rfc[i] to the molcule's flux.
         here we include the sign with the stoichiometric_coef, so for
	 reactants the rfc contribution is subtracted, and for products it
	 is added.
  */
  for (i=0;i<num_rxns;i++) {
    pt = 1.0;
    rt = 1.0;
    tr = 1.0;
    tp = 1.0;
    multiplier = 1.0;
    for (j=rxn_ptrs[i];j<rxn_ptrs[i+1];j++) {
      mi = molecule_indices[j];
      molecule = (struct molecule_struct *)&molecules[mi];
      ci = molecule->c_index;
      compartment = (struct compartment_struct *)&compartments[ci];
      recip_volume       = compartment->recip_volume;
      volume             = compartment->volume;
      volume_avo         = volume * avogadro;
      recip_volume_avo   = recip_volume * recip_avogadro;
      /*
      */
      klim = coefficients[j];
      conc_mi = concs[mi];
      if (klim < 0.0) {
	klim = 0.0 - klim;
	rt = rt * conc_to_pow(conc_mi,klim,factorial);
	multiplier = multiplier * conc_to_pow(volume_avo,klim,factorial);
	/*
	for (k=0;k<(-klim);k++) {
	  rt = rt * conc_mi;
	  multiplier *= volume;
	  tr = tr * (count_mi - klim);
	}
	*/
      } else {
	if (klim > 0.0) {
	  pt = pt * conc_to_pow(conc_mi,klim,factorial);
	  multiplier = multiplier * conc_to_pow(recip_volume_avo,klim,
						factorial);
	  /*
	  for (k=0;k<klim;k++) {
	    pt = pt * conc_mi;
	    multiplier *= recip_volume;
	    tp = tp * (count_mi + klim);
	  }
	  */
	}
      }
    } /* end for (j...) */
    keq_adj = multiplier;
    rkeq_adj = 1.0;
    /*
    sum_coeff = coeff_sum[i];
    if (sum_coeff > 0) {
      for (k=0;k < sum_coeff; k++) {
	keq_adj *= recip_avogadro;
      }
    } else {
      if (sum_coeff < 0) {
	sum_coeff = - sum_coeff;
	for (k=0;k < sum_coeff; k++) {
	  keq_adj *= avogadro;
	}
      }
    }
    */
    if (keq_adj != 0.0) {
      rkeq_adj = 1.0/keq_adj;
    }
    /*
      NB. tp and tr will always be > 0 as |klim| > 0 and concs_mi >= 0;
      but now the reaction contribution is in counts per time, if we
      want it in moles/liter/time we need to divide by volume and
      avogadro's number - wonder if we really need these two multiplies??
    rfc[i] = (ke[i] * (rt/tp)) - (rke[i] * (pt/tr)) * recip_volume * recip_avogadro;
    */
    /*
      Save likelihoods for printing.
    */
    if (pt != 0) {
      forward_lklhd[i] = ke[i] * keq_adj * (rt/pt);
    } else {
      if (rt != 0) {
	forward_lklhd[i] = 1.0;
      } else {
	forward_lklhd[i] = 0.0;
      }
    }
    if (rt != 0) {
      reverse_lklhd[i] = rke[i] * rkeq_adj * (pt/rt);
    } else {
      if (pt != 0) {
	reverse_lklhd[i] = 1.0;
      } else {
	reverse_lklhd[i] = 0.0;
      }
    }
    /*
    rfc[i] = (ke[i] * (rt/tp)) - (rke[i] * (pt/tr));
    NB if use_activities is not set activities[i] will be 1.0 for all i.
    */
    rfc[i] = (forward_lklhd[i] - reverse_lklhd[i]) * activities[i];
  } /* end for i */
  if (success) {
    molecule = molecules;
    for (i=0;i<num_species;i++) {
      fluxi = 0.0;
      if (molecule->variable == 1) {
	for (j=molecules_ptrs[i];j<molecules_ptrs[i+1];j++) {
	  rxn = rxn_indices[j];
	  coeff = coefficients[j];
	  if (coeff != 0.0) {
	    fluxi += (rfc[rxn]*coeff);
	  }
	} /* end for(j...) */
	flux[i] = flux_scaling * fluxi;
      } else {
	flux[i] = 0.0;
      }
      molecule += 1; /* Caution address arithmetic here. */
    } /* end for (i...) */
  } /* end if success */
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"Mol_index\t   conc   \t    flux\n");
    for (i=0;i<num_species;i++) {
      fprintf(lfp,"%d\t%le\t%le\n",
	      i,concs[i],flux[i]);
    }
    fflush(lfp);
  }
#endif
  return (success);
}
