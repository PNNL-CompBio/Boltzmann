#include "boltzmann_structs.h"
#include "boltzmann_cvodes_headers.h"
#include "cvodes_params_struct.h"
#include "conc_to_pow.h"
#include "update_regulations.h"
#include "lr7_gradient.h"

int lr7_gradient(struct state_struct *state, 
		 double *concs,
		 double *flux, 
		 int choice) {
  /*
    Compute approximations to concentration changes wrt time, 
    Based on thermodynamics formulation for concentration rates of change.

    Get reference from Bill Cannon
    Called by: gradient
    Calls:     update_regulations, conc_to_pow, fprintf, fflush

                                TMF
    state                       *SI   Boltzmant state structure.
                                      uses number_reactions,
				           unique_moleules,
					   sorted_molecules,
					   sorted_compartments,
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
  double  *rfc;
  double  *ke;
  double  *rke;
  double  *coefficients;
  double  *rcoefficients;
  double  *coeff_sum;
  double  *kq;
  double  *kqi;

  int64_t *molecules_ptrs;
  int64_t *rxn_indices;
  int64_t *rxn_ptrs;
  int64_t *molecule_indices;

  double  flux_scaling;
  double  pt;
  double  rt;
  double  tr;
  double  tp;
  double  conc_mi;
  /*
  double  thermo_adj;
  */
  double  volume;
  double  recip_volume;
  double  avogadro;
  double  recip_avogadro;
  double  fluxi;
  double  keq_adj;
  double  rkeq_adj;
  double  multiplier;
  double  klim;
  double  factorial;
  double  volume_avo;
  double  recip_volume_avo;
  double  conc_mi_plus;
  /*
  double  sum_coeff;
  */

  int num_species;
  int num_rxns;

  int rxn;
  int success;

  int i;
  int j;

  int mi;
  int ci;

  int compute_sensitivities;
  int ode_solver_choice;

  int use_regulation;
  int count_or_conc;

  FILE *lfp;
  FILE *efp;
  /*
#define DBG 1
  */
  /*
    Check that base_rxn is in range.
  */
  success = 1;
  num_rxns         = state->number_reactions;
  num_species      = state->nunique_molecules;
  molecules        = state->sorted_molecules;
  activities       = state->activities;
  compartments     = state->sorted_compartments;
  molecules_matrix = state->molecules_matrix;
  coeff_sum        = state->coeff_sum;
  molecules_ptrs   = molecules_matrix->molecules_ptrs;
  rxn_indices      = molecules_matrix->reaction_indices;
  coefficients     = molecules_matrix->coefficients;
  rxn_matrix       = state->reactions_matrix;
  rxn_ptrs         = rxn_matrix->rxn_ptrs;
  molecule_indices = rxn_matrix->molecules_indices;
  rcoefficients    = rxn_matrix->coefficients;
  ke               = state->ke;
  rke              = state->rke;
  rfc              = state->product_term;
  avogadro         = state->avogadro;
  recip_avogadro   = state->recip_avogadro;
  use_regulation   = state->use_regulation;
  kq               = state->ode_kq;
  kqi              = state->ode_kqi;
  factorial      = 0.0;
  /*
    If we are using cvodes and computing sensitivites the 
    call may be made with perturbed equilibrium constants (the sensitivity
    parameters), so take them from the cvodes_params vector.
  */
  ode_solver_choice = state->ode_solver_choice;
  compute_sensitivities = state->compute_sensitivities;
  if ((ode_solver_choice == 1) && compute_sensitivities) {
    cvodes_params = state->cvodes_params;
    ke = cvodes_params->p;
    rke = cvodes_params->rp;
    for (i=0;i<num_rxns;i++) {
      rke[i] = 1.0/ke[i];
    }
  }
  /*
  flux_scaling     = compute_flux_scaling(state,concs);
  */
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
            k_r * product of products*stoichiometric_coef,

	 where k_f = k_eq/ thermo_product of products.

	 and   k_r  = k_eq^(-1)/thermo_product of reactants,

	 and thermo_product = 
	 product( species_conc + |stoichiometric_coef|/volume)^|stoichiometric coef|.

	 Then if a molecule is a  in rxn i,

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
    /*
    sum_coeff = 0.0;
    */
    for (j=rxn_ptrs[i];j<rxn_ptrs[i+1];j++) {
      mi = molecule_indices[j];
      molecule = (struct molecule_struct *)&molecules[mi];
      ci = molecules->c_index;
      compartment = (struct compartment_struct *)&compartments[ci];
      recip_volume = compartment->recip_volume;
      volume       = compartment->volume;
      volume_avo   = volume * avogadro;
      recip_volume_avo = recip_volume * recip_avogadro;
      klim = rcoefficients[j];
      /*
      sum_coeff += klim;
      thermo_adj = abs(klim) * recip_volume_avo;
      */
      conc_mi = concs[mi];
      if (klim < 0.0) {
	klim = 0.0 - klim;
	conc_mi_plus = conc_mi + (klim * recip_volume_avo);
	rt = rt * conc_to_pow(conc_mi,klim,factorial);
	tr = tr * conc_to_pow(conc_mi_plus,klim,factorial);
	multiplier = multiplier * conc_to_pow(volume_avo,klim,factorial);
	/*
	for (k=0;k<(-klim);k++) {
	  rt = rt * conc_mi;
	  tr = tr * (conc_mi + thermo_adj);
	}
	*/
      } else {
	if (klim > 0.0) {
	  conc_mi_plus = conc_mi + (klim * recip_volume_avo);
	  pt = pt * conc_to_pow(conc_mi,klim,factorial);
	  tp = tp * conc_to_pow(conc_mi_plus,klim,factorial);
	  multiplier = multiplier * conc_to_pow(recip_volume_avo,klim,
						factorial);
	  /*
	  for (k=0;k<klim;k++) {
	    pt = pt * conc_mi;
	    tp = tp * (conc_mi + thermo_adj);
	  }
	  */
	}
      }
    } /*  end for (j... ) */
    keq_adj = multiplier;
    rkeq_adj = 1.0/multiplier;
    /*
    keq_adj = 1.0;
    rkeq_adj = 1.0;
    multiplier = 1.0;
    sum_coeff = coeff_sum[i];
    if (sum_coeff > 0) {
      multiplier = recip_volume * recip_avogadro;
    } else {
      if (sum_coeff < 0) {
	multiplier = volume * avogadro;
	sum_coeff = - sum_coeff;
      } 
    }
    for (k=0;k < sum_coeff; k++) {
      keq_adj *= multiplier;
    }
    rkeq_adj = 1.0/keq_adj;
    */
    /*
    fprintf(stderr,"for reaction %d keq_adj = %le", i,keq_adj);
    */
    /*
      NB. tp and tr will always be > 0 as thermo_adj > 0 and concs_mi >= 0;
      thermo_adj and 1/coefficients[j] could be precomputed, and stored
      as fields in the reaction and molecules matrices respectively.
      NB if use_activities is not set activities[i] will be 1.0 for all i.
    */
    kq[i] = keq_adj * ke[i] * (rt/tp);
    kqi[i] = rkeq_adj * rke[i] * pt/tr;
    rfc[i] = ((ke[i] * keq_adj * (rt/tp)) - (rke[i] * rkeq_adj * (pt/tr))) * activities[i];
    /*
    fprintf(stderr,"Rxn %d: KQ = %e K = %e  Q^-1 = %e  keq_adj = %e  activities = %e\n", i, rfc[i], ke[i], rt/tp, keq_adj, activities[i]);
    */
  }
  if (success) {
    molecule = molecules;
    for (i=0;i<num_species;i++) {
      fluxi = 0.0;
      if (molecule->variable == 1) {
	for (j=molecules_ptrs[i];j<molecules_ptrs[i+1];j++) {
	  rxn = rxn_indices[j];
	  if (coefficients[j] != 0.0) {
	    fluxi += (rfc[rxn]*((double)coefficients[j]));
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
