#include "boltzmann_structs.h"
#include "get_counts.h"
#include "lr8_approximate_ys0.h"
int lr8_approximate_ys0(struct state_struct *state, double *ys0v,
			double *concs) {
  /*
    Approximate df/dke
    Called by: approximate_ys0
  */
  struct  molecules_matrix_struct *molecules_matrix;
  struct  reactions_matrix_struct *rxn_matrix;
  double *activities;
  double *ys0vi;
  double *ke;
  double *rke;
  double  pt;
  double  rt;
  double  tr;
  double  tp;
  double  *counts;
  double  *conc_to_count;
  double  drfc;
  int64_t *rxn_ptrs;
  int64_t *molecules_indices;
  int64_t *rcoefficients;
  int64_t *molecules_ptrs;
  int64_t *rxn_indices;
  int64_t *mcoefficients;
  int ny;
  int ns;

  int i;
  int j;

  int success;
  int mi;

  int klim;
  int k;

  int count_mi;
  int padi;

  success           = 1;
  ke  		    = state->ke;
  rke 		    = state->rke;
  ny  		    = state->nunique_molecules;
  ns  		    = state->number_reactions;
  counts            = state->ode_counts;
  conc_to_count     = state->conc_to_count;
  activities        = state->activities;
  rxn_matrix        = state->reactions_matrix;
  rxn_ptrs          = rxn_matrix->rxn_ptrs;
  molecules_indices = rxn_matrix->molecules_indices;
  rcoefficients     = rxn_matrix->coefficients;
  molecules_matrix  = state->molecules_matrix;
  molecules_ptrs    = molecules_matrix->molecules_ptrs;
  rxn_indices       = molecules_matrix->reaction_indices;
  mcoefficients     = molecules_matrix->coefficients;
  get_counts(ny,concs,conc_to_count,counts);
  ys0vi = ys0v;
  for (i=0;i<ns;i++) {
    pt = 1.0;
    rt = 1.0;
    tr = 1.0;
    tp = 1.0;
    for (j=rxn_ptrs[i];j<rxn_ptrs[i+1];j++) {
      mi = molecules_indices[j];
      /*
      molecule = (struct molecule_struct *)&molecules[mi];
      ci = molecule->c_index;
      compartment = (struct compartment_struct *)&compartments[ci];
      recip_volume       = compartment->recip_volume;
      */
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
    drfc = ((rt/tp) + (rke[i]*rke[i]*(pt/tr))) * activities[i];
    for (j=0;j<ny;j++) {
      ys0vi[j] = 0.0;
    }
    for (j=rxn_ptrs[i];j<rxn_ptrs[i+1];j++) {
      mi = molecules_indices[j];
      ys0vi[mi] = rcoefficients[j]*drfc;
    }
    ys0vi += ny; /* Caution address arithmetic */
  } /* end for i */
  return(success);
}
