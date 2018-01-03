#include "boltzmann_structs.h"
#include "fill_flux_pieces.h"
int fill_flux_pieces(struct state_struct *state, 
		     struct molecules_matrix_struct *molecules_matrix,
		     int base_reaction,
		     int fill_jacobi) {
  /*
    Compute reactant_terms, product_terms, p_over_r, r_over_p,
    flux_vector and flux_jacobian outputs from inputs,
    ke, fwrd_lklhd, rvrs_lklhd, concs, the reactions_matrix and the
    molecules_matrix.
    Called by: deq_run
  */
  struct rxn_matrix_struct *reaction_matrix;
  struct molecule_struct *molecules;
  struct molecule_struct *molecule;
  struct compartment_struct *compartments;
  struct compartment_struct *compartment;
  double *ke;
  double *fwrd_lklhd;
  double *rvrs_lklhd;
  double *concs;
  double *reactant_term;
  double *product_term;
  double *p_over_r;
  double *r_over_p;
  double *fluxes;
  double *flux_jacobian;
  double *jacobian_row;
  double p;
  double r;
  double q;
  double rq;
  double lklhd_b;
  double lklhd_rb;
  double recip_lklhd_b;
  double conc;
  double conc_recip;
  double flux;
  double fluxfluxf;
  double fluxfluxr;
  int64_t *rxn_ptrs;
  int64_t *mcoefficients;
  int64_t *m_index;
  int64_t *mol_ptrs;
  int64_t *r_index;
  int64_t *rcoefficients;
  int    num_species;
  int    num_rxns;

  int    success;
  int    i;

  int    j;
  int    k;

  int    mindex;
  int    alpha;

  int    rindex;
  int    cindex;

  int    beta;
  int    padi;

  FILE *lfp;
  FILE *efp;
  success          = 1;
  lfp              = state->lfp;
  ke               = state->ke;
  fwrd_lklhd       = state->forward_rxn_likelihood;
  rvrs_lklhd       = state->reverse_rxn_likelihood;
  concs            = state->concs;
  reaction_matrix  = state->reactions_matrix;
  reactant_term    = state->reactant_term;
  product_term     = state->product_term;
  p_over_r         = state->p_over_r;
  r_over_p         = state->r_over_p;
  fluxes           = state->flux_vector;
  flux_jacobian    = state->flux_jacobian;
  num_rxns         = (int)state->number_reactions;
  num_species      = (int)state->nunique_molecules;
  rxn_ptrs         = reaction_matrix->rxn_ptrs;
  m_index          = reaction_matrix->molecules_indices;
  mcoefficients    = reaction_matrix->coefficients;
  mol_ptrs         = molecules_matrix->molecules_ptrs;
  r_index          = molecules_matrix->rxn_indices;
  rcoefficients    = molecules_matrix->coefficients;
  molecules        = state->sorted_molecules;
  compartments     = state->sorted_cmpts;
  if (base_reaction >= 0) {
    if (base_reaction < num_rxns) {
      lklhd_b = fwrd_lklhd[base_reaction];
      lklhd_rb = rvrs_lklhd[base_reaction];
      if (lklhd_b > 0) {
	recip_lklhd_b = 1.0/lklhd_b;
      } else {
	success = 0;
      }
    } else {
      if (base_reaction < num_rxns+num_rxns) {
	lklhd_b = rvrs_lklhd[base_reaction-num_rxns];
	lklhd_rb = fwrd_lklhd[base_reaction-num_rxns];
	if (lklhd_b > 0) {
	  recip_lklhd_b = 1.0/lklhd_b;
	} else {
	  success = 0;
	}
      } else {
	success = 0;
      }
    }
  } else {
    success = 0;
  }
  if (success == 0) {
    if (lfp) {
      fprintf(lfp,"fill_flux_pieces: error invalid base_reaction number was %d\n",
	      base_reaction);
      fflush(lfp);
    }
  }
  if (success) {
    /*
      Loop over each reaction computing reactant_term and product_term,
      p_over_r, r_over_p,
    */
    for (i=0;i<num_rxns;i++) {
      r = 1.0;
      p = 1.0;
      for (j=rxn_ptrs[i];j<rxn_ptrs[i+1];j++) {
	mindex = m_index[j];
	alpha    = (int)mcoefficients[j];
	conc     = concs[mindex];
	if (alpha < 0) {
	  for (k=0;k<(-alpha);k++) {
	    r = r * conc;
	  }
	} else {
	  for (k=0;k<alpha;k++) {
	    p = p * conc;
	  }
	}
      }
      reactant_term[i] = r;
      product_term[i]  = p;
      p_over_r[i]       = 0.0;
      r_over_p[i]       = 0.0;
      if (r > 0) {
	p_over_r[i]     = p/r;
      }
      if (p > 0) {
	r_over_p[i]     = r/p;
      }
    }
  }
  /*
    Now compute fluxes for each molecule summing over reactions.
  */
  if (success) {
    for (i=0;i<num_species;i++) {
      flux = 0.0;
      for (j=mol_ptrs[i];j<mol_ptrs[i+1];j++) {
	rindex = r_index[j];
	alpha    = (int)rcoefficients[j];
	r        = reactant_term[rindex];
	p        = product_term[rindex];
	q        = p_over_r[rindex];
	rq       = r_over_p[rindex];
	if (alpha < 0) {
	  /*
	    Species i is a reactant in reaction.
	  */
	  if (r > 0.0) {
	    flux = flux - fwrd_lklhd[rindex]*recip_lklhd_b*(1.0 - (q/ke[rindex]));
	  } else {
	    if (p > 0.0) {
	      flux = flux + rvrs_lklhd[rindex]*recip_lklhd_b*(1.0 - (rq * ke[rindex]));
	    }
	  }
	} else {
	  /*
	    Species is a product in reaction.
	  */
	  if (r > 0.0) {
	    flux = flux + fwrd_lklhd[rindex]*recip_lklhd_b*(1.0 - (q/ke[rindex]));
	  } else {
	    if (p > 0.0) {
	      flux = flux - rvrs_lklhd[rindex]*recip_lklhd_b*(1.0 - (rq * ke[rindex]));
	    }
	  }
	}
      } /* end for (j...) */
      fluxes[i] = flux;
    } /* end for (i...) */
  }
  if (success && fill_jacobi) {
    for (i=0;i<num_species*num_species;i++) {
      flux_jacobian[i] = 0.0;
    }
    jacobian_row     = flux_jacobian;
    for (i=0;i<num_species;i++) {
      for (j=mol_ptrs[i];j<mol_ptrs[i+1];j++) {
	rindex = r_index[j];
	alpha    = (int)rcoefficients[j];
	r        = reactant_term[rindex];
	p        = product_term[rindex];
	q        = p_over_r[rindex];
	rq       = r_over_p[rindex];
	if (alpha < 0) {
	  if (r > 0) {
	    fluxfluxf = - fwrd_lklhd[rindex] * recip_lklhd_b;
	  } else {
	    fluxfluxf = 0.0;
	  }
	  if (p > 0) {
	    fluxfluxr =   rvrs_lklhd[rindex] * recip_lklhd_b;
	  } else {
	    fluxfluxr = 0.0;
	  }
	} else {
	  if (r > 0) {
	    fluxfluxf = fwrd_lklhd[rindex] * recip_lklhd_b;
	  } else {
	    fluxfluxf = 0.0;
	  }
	  if (p > 0) {
	    fluxfluxr =  -rvrs_lklhd[rindex] * recip_lklhd_b;
	  } else {
	    fluxfluxr = 0.0;
	  }
	}
	for (k=rxn_ptrs[rindex];k<rxn_ptrs[rindex+1];k++) {
	  mindex = m_index[k];
	  conc   = concs[mindex];
	  beta   = mcoefficients[k];
	  if (conc > 0.0) {
	    conc_recip = 1.0/conc;
	    if (beta < 0) {
	      jacobian_row[mindex] = fluxfluxf * conc_recip;
	    } else {
	      jacobian_row[mindex] = fluxfluxr * conc_recip;
	    }
	  } /* end if conc > 0 */
	} /* end for (k...) */
      } /* end for (j...) */
      jacobian_row = jacobian_row + num_species; /* Caution address arithmetic */
    } /* end for (i...) */
  }
  return(success);
}

