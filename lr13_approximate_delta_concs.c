#include "boltzmann_structs.h"
#include "get_counts.h"
#include "update_regulations.h"
#include "lr13_approximate_delta_concs.h"

int lr13_approximate_delta_concs(struct state_struct *state, 
				double *concs,
				double *flux, 
				int choice) {
  /*
    Compute approximations to concentration changes wrt time, 
    based on entropy?

    d n_i/dt = 1/log(n_i/exp(-molecule_dg0tfs[i]/RT)) sum K_j *  Q_j^-1
                                                       j

    n_i  = count for species i.

    Q_j^-1 = reactants/products

    in reactions j where molecule i is involved.
      
    K_j = keq[j]

    Get reference from Bill Cannon
    Called by: approximate_delta_concs
    Calls:     get_counts,
               update_regulations, fabs, log, exp

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
  struct  molecule_struct *molecules;
  struct  molecule_struct *molecule;
  /*
  struct  compartment_struct *compartments;
  struct  compartment_struct *compartment;
  */
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
  double  *molecule_dg0tfs;
  double  flux_scaling;
  double  pt;
  double  rt;
  double  tr;
  double  tp;
  double  dgi;
  double  qii;
  double  m_r_rt;
  /*
  double  conc_mi;
  double  thermo_adj;
  double  recip_volume;
  double  recip_avogadro;
  */
  double  fluxi;
  double  coef;
  double  count_mi;
  int64_t *molecules_ptrs;
  int64_t *rxn_indices;
  int64_t *coefficients;
  int64_t *rxn_ptrs;
  int64_t *molecule_indices;
  int64_t *rcoefficients;
  int num_species;
  int num_rxns;

  int rxn;
  int success;

  int i;
  int j;

  int mi;
  int padi;

  int k;
  int klim;

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
  success          = 1;
  num_rxns         = state->number_reactions;
  num_species      = state->nunique_molecules;
  molecules        = state->sorted_molecules;
  /*
  compartments     = state->sorted_compartments;
  */
  activities       = state->activities;
  molecule_dg0tfs  = state->molecule_dg0tfs;
  forward_lklhd    = state->ode_forward_lklhds;
  reverse_lklhd    = state->ode_reverse_lklhds;
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
  m_r_rt           = state->m_r_rt;
  /*
  recip_avogadro   = state->recip_avogadro;
  */
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

  */
  for (i=0;i<num_rxns;i++) {
    pt = 1.0;
    rt = 1.0;
    tr = 1.0;
    tp = 1.0;
    for (j=rxn_ptrs[i];j<rxn_ptrs[i+1];j++) {
      mi = molecule_indices[j];
      klim = rcoefficients[j];
      count_mi = counts[mi];
      if (klim < 0) {
	for (k=0;k<(-klim);k++) {
	  rt = rt * count_mi;
	  tr = tr * (count_mi - klim);
	}
      } else {
	if (klim > 0) {
	  for (k=0;k<klim;k++) {
	    pt = pt * count_mi;
	    tp = tp * (count_mi + klim);
	  }
	}
      }
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
    forward_lklhd[i] = ke[i] * (rt/tp);
    reverse_lklhd[i] = rke[i] * (pt/tr);
    /*
    rfc[i] = (ke[i] * (rt/tp)) - (rke[i] * (pt/tr));
    NB if use_activities is not set activities[i] will be 1.0 for all i.
    */
    /*
    if (pt != 0.0) {
      qii = rt/pt;
      rfc[i] = log(ke[i]*qii);    
    } else {
      rfc[i] = 0.0;
      if (lfp) {
	fprintf(lfp,"lr13_approximate_delta_concs pt for reacion %d was 0.0. Setting rfc[%d] to 0.0\n",i,i);
	fflush(lfp);
      }
    }
    */
    qii = rt/tp;
    rfc[i] = log(ke[i]*qii);    
  }
  if (success) {
    molecule = molecules;
    for (i=0;i<num_species;i++) {
      fluxi = 0.0;
      count_mi = counts[i];
      if (molecule->variable == 1) {
	if (count_mi <= 0.0) {
	  count_mi = .00001;
	  if (lfp) {
	    fprintf(lfp,"lr13_approximate_delta_concs: 0 count for molecule %d was set to .00001\n",i);
	    fflush(lfp);
	  }
	}
	dgi = molecule_dg0tfs[i];
	flux_scaling = 1.0/log(count_mi/exp(dgi*m_r_rt));
	for (j=molecules_ptrs[i];j<molecules_ptrs[i+1];j++) {
	  rxn = rxn_indices[j];
	  coef = (double)coefficients[j];
	  if (coef != 0.0) {
	    fluxi += rfc[rxn]/coef;
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
