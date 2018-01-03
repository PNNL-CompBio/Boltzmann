#include "boltzmann_structs.h"
#include "get_counts.h"
#include "update_regulations.h"
#include "vec_set_constant.h"
#include "crs_column_sort_rows.h"
#include "lr8_approximate_jacobian.h"

int lr8_approximate_jacobian(struct state_struct *state, 
			 double *concs,
			 double *delta_concs,
			 double t,
			 int choice) {
			 
  /*
    Compute approximations to the jacobian of theconcentration changes. wrt time,   Based on thermodynamics formulation for concentraion rate
    changes     using counts to compute tr, tp, pt, rt instead of concs.
    Based on lr8_approximate delta concs.
    drfc could be taken from state and is a number_reactions * nunique_molecules
    matrix.
    dfdy is the jacbian output and is nunqiue_molecules * nunique_molecules,
    may eventually want to put this out in compressed row storage format,
    as it is a sparse matrix:  ia_dfdy, ja_dfdy, dfdy
    Also might want to make those fields in state instead of passing in
    as arguments.


    Called by: approximate_jacobian
    Calls:     get_counts,
               update_regulations,
	       vec_set_constant,
	       crs_column_sort_rows

                                TMF
    state                       *SI   Boltzmant state structure.
                                      uses number_reactions,
				           unique_moleules,
					   sorted_molecules,
                                           molecules_matrix,
					   ke, rke,
					   product_term as scratch.
					   and lfp,
					   drfc[num_molecules]
					   dfdy_row[ny]
					   dfdy_a[<= ny*ny]
					   dfdy_at<= ny*ny]
					   drfc[num_molecules]
					   dfdy_ja[<-ny*ny]
					   dfdy_jat[<-ny*ny]
					   dfdy_ia[ny+1]
					   dfdy_iat[ny+1]

    concs			D1I   molecule concentrations vector of length 
                                      nunique_moleucles

    drfc                        D1O   vector of length number_reactions * 
                                      unique_molecules
                                      of partial derivatives of 
				      concentration change with respect
				      to other concentrations.
				      Set by this routine.

    choice                      IOI   Not used by this routine.

  */
  struct  molecule_struct *molecules;
  struct  molecule_struct *molecule;
  struct  molecules_matrix_struct *molecules_matrix;
  struct  reactions_matrix_struct *rxn_matrix;
  struct  cvodes_params_struct *cvodes_params;
  double  *activities;
  double  *forward_lklhd;
  double  *reverse_lklhd;
  double  *rfc;
  double  *ke;
  double  *rke;
  double  *counts;
  double  *conc_to_count;
  double  *drfc;
  double  *dfdy_a;
  double  *dfdy_at;
  double  *dfdy_row;
  double  *recip_coeffs;
  double  flux_scaling;
  double  pt;
  double  rt;
  double  tr;
  double  tp;
  double  flklhd;
  double  rlklhd;
  double  recip_coefficient;
  /*
  double  conc_mi;
  double  thermo_adj;
  double  recip_volume;
  double  recip_avogadro;
  */
  double  fluxi;
  double  count_mi;
  double  activityi;
  double  dzero;
  int64_t *molecules_ptrs;
  int64_t *rxn_indices;
  int64_t *coefficients;
  int64_t *rxn_ptrs;
  int64_t *molecule_indices;
  int64_t *rcoefficients;
  int     *dfdy_ja;
  int     *dfdy_ia;
  int     *dfdy_jat;
  int     *dfdy_iat;
  int     *column_mask;
  int ny;
  int num_rxns;

  int rxn;
  int success;

  int i;
  int j;

  int mi;
  int dfdy_pos;

  int k;
  int coef;

  int use_regulation;
  int count_or_conc;

  
  int mk;
  int mj;

  FILE *lfp;
  FILE *efp;
  /*
#define DBG 1
  */
  success          = 1;
  dzero            = 0.0;
  num_rxns         = state->number_reactions;
  ny      = state->nunique_molecules;
  molecules        = state->sorted_molecules;
  activities       = state->activities;
  forward_lklhd    = state->ode_forward_lklhds;
  reverse_lklhd    = state->ode_reverse_lklhds;
  cvodes_params    = state->cvodes_params;
  molecules_matrix = state->molecules_matrix;
  molecules_ptrs   = molecules_matrix->molecules_ptrs;
  rxn_indices      = molecules_matrix->reaction_indices;
  coefficients     = molecules_matrix->coefficients;
  recip_coeffs     = molecules_matrix->recip_coeffs;
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
  /*
    The following vectors are allocated in boltzmann_cvodes 
  */
  drfc             = cvodes_params->drfc;    
  dfdy_a           = cvodes_params->dfdy_a;  
  dfdy_ia          = cvodes_params->dfdy_ia; 
  dfdy_ja          = cvodes_params->dfdy_ja; 
  dfdy_at          = cvodes_params->dfdy_at; 
  dfdy_iat         = cvodes_params->dfdy_iat; 
  dfdy_jat         = cvodes_params->dfdy_jat; 
  dfdy_row         = cvodes_params->prec_row;
  column_mask      = cvodes_params->column_mask;
  /*
  recip_avogadro   = state->recip_avogadro;
  */
  /*
  flux_scaling     = compute_flux_scaling(state,concs);
  */
  get_counts(ny,concs,conc_to_count,counts);
  flux_scaling     = 1.0;
  lfp      = state->lfp;
  /*
    As per discusion with Bill Cannon, we want to update the activities
    if regulation is in play. So do that here.
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
    activityi = activities[i];
    for (j=rxn_ptrs[i];j<rxn_ptrs[i+1];j++) {
      mi = molecule_indices[j];
      coef = (int)rcoefficients[j];
      count_mi = counts[mi];
      if (coef < 0) {
	for (k=0;k<(-coef);k++) {
	  rt = rt * count_mi;
	  tr = tr * (count_mi - coef);
	}
      } else {
	if (coef > 0) {
	  for (k=0;k<coef;k++) {
	    pt = pt * count_mi;
	    tp = tp * (count_mi + coef);
	  }
	}
      }
    } /*end for (j...) */
    /*
      NB. tp and tr will always be > 0 as |coef| > 0 and concs_mi >= 0;
      but now the reaction contribution is in counts per time, if we
      want it in moles/liter/time we need to divide by volume and
      avogadro's number - wonder if we really need these two multiplies??
    rfc[i] = (ke[i] * (rt/tp)) - (rke[i] * (pt/tr)) * recip_volume * recip_avogadro;
    */
    flklhd = ke[i] * (rt/tp);
    rlklhd = rke[i] * (pt/tr);
    /*
    rfc[i] = (ke[i] * (rt/tp)) - (rke[i] * (pt/tr)) * activites[i];
    NB if use_activities is not set activities[i] will be 1.0 for all i.
    */
    for (j=rxn_ptrs[i];j<rxn_ptrs[i+1];j++) {
      mi = molecule_indices[j];
      molecule = (struct molecule_struct *)&molecules[mi];
      if (molecule->variable) {
	/*
	  d/dy_mi for reaction i
	*/
	coef = (int)rcoefficients[j];
	count_mi = counts[mi];
	if (coef < 0) {
	  if (count_mi > 0) {
	    drfc[j] =  0.0 -coef * ((flklhd/count_mi) + (rlklhd/(count_mi - coef))) * activityi;	
	  } else {
	    drfc[j] = rlklhd * activityi;
	  }
	} else {
	  if (coef > 0) {
	    if (count_mi > 0) {
	      drfc[j] = 0.0 - coef *((flklhd/(count_mi + coef)) + (rlklhd/count_mi)) * activityi;
	    } else {
	      drfc[j] = - flklhd * activityi;
	    }
	  } else {
	    drfc[j] = 0.0;
	  }
	}
      } else {
	drfc[j] = 0.0;
      } /* end if (molecule->variable) */
    } /* end for j */
  } /* end for i */
  molecule = molecules;
  dfdy_pos = 0;
  for (i=0;i<ny;i++) {
    column_mask[i] = 0;
  }
  dfdy_ia[0]  = 0;
  for (i=0;i<ny;i++) {
    fluxi = 0.0;
    vec_set_constant(ny,dfdy_row,dzero);
    /*
      Ensure that we store a diagonal element for each row.
    */
    column_mask[i] = 1;
    dfdy_ja[dfdy_pos] = i;
    dfdy_pos += 1;
    if (molecule->variable == 1) {
  	for (j=molecules_ptrs[i];j<molecules_ptrs[i+1];j++) {
  	  rxn = rxn_indices[j];
  	  if (coefficients[j]  != 0) {
  	    recip_coefficient = recip_coeffs[j];
  	    /*
  	      Add rxn row of drfc * recip_coefficient to row i of dfdy
  	    */
  	    for (k=rxn_ptrs[rxn];k<rxn_ptrs[rxn+1];k++) {
  	      mk = molecule_indices[k];
  	      if (column_mask[mk] == 0) {
  		column_mask[mk] = 1;
  		dfdy_ja[dfdy_pos] = mk;
  		dfdy_pos += 1;
  	      }
  	      dfdy_row[mk] += drfc[k] * recip_coefficient;
  	    }
  	  }
  	}
  	/*
  	  The column numbers of nonzero elements in row i of
  	  dfdy are now in dfdy_ja[dfdy_ia[i]:dfdy_pos-1]
  	  Ultimately we want to sort them to so as 
  	  to make preconditioning computation and mvp computation 
  	  more effiecient, but that is more efficiently done
  	  with a double transpose algorithm.

  	  Extract the sparse dfdy row from the dfd_row vector,
  	  reseting it and the column_mask vector as we go.
  	*/
  	for (j=dfdy_ia[i];j<dfdy_pos;j++) {
  	  mj = dfdy_ja[j];
  	  dfdy_a[j] = dfdy_row[mj];
  	  dfdy_row[mj] = 0.0;
  	  column_mask[mj] = 0;
  	}
    } /* end if (molecule->variable) */
    dfdy_ia[i+1] = dfdy_pos;
  } /* end for (i...) */
  /*
    Now we want rows of dfdy in dfdy_a to be sorted by column number,
    This is best accomplished with a double transpose algorithm
    requiring work space = sizeof (dfdy_a) + sizeof(dfdy_ja) + sizeof (dfdy_ia).
  */
  crs_column_sort_rows(ny, ny, dfdy_a, dfdy_ia, dfdy_ja, dfdy_at, dfdy_iat, 
		       dfdy_jat);

  return (success);
}
