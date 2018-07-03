#include "boltzmann_structs.h"
#include "vec_set_constant.h"
#include "conc_to_pow.h"
#include "compute_dfdke_dfdmu0.h"
int compute_dfdke_dfdmu0(struct state_struct *state, double *concs) {
  /*
    Compute the partials of the conentration deriviatives with respect to
    the equilibrium constants and the chemical potential creating the
    .dfdke and .dfdmu0 files at the end of ode simulation.
    scratch needs to be a vectro of length at least 3*nunique_molecules

    Called by ode_solver.c
    Calls     vec_set_constant, conc_to_pow, fopen, fprintf, fclose
  */
  /*
    Open the .dfdke and dfdmu0 files.
    Compute the pieces of the partials, the forward and reverse pieces so
    that if we desire we may add them stably for the partials wrt the
    chemical potentials.  As these matrices of partials are sparse,
    we only print the nonzero paritials.
    
  */
  struct  molecules_matrix_struct *molecules_matrix;
  struct  reactions_matrix_struct *rxn_matrix;
  double  *activities;
  double  *forward_piece;
  double  *reverse_piece;
  double  *mu0_partials;
  double  *ke;
  double  *rke;
  double  *counts;
  double  *conc_to_count;
  double  *dfdke_dfdmu0_work;
  int64_t *molecules_ptrs;
  int64_t *molecule_indices;
  int64_t *rxn_ptrs;
  int64_t *rxn_indices;
  double  *rcoefficients;
  double  *mcoefficients;
  char    *dfdke_file;
  char    *dfdmu0_file;

  double  count_mi;
  double  pt;
  double  rt;
  double  tp;
  double  tr;
  double  dzero;
  double  multiplier;
  double  m_r_rt;
  double  gamma_im;
  double  gamma_ik;
  double  dfdke;
  double  count_plus;
  double  telescoping;

  int num_species;
  int num_rxns;

  int i;
  int j;

  int success;
  int mi;

  int rj;
  int delta_concs_choice;

  int k;
  int mk;



  FILE *dfdke_fp;
  FILE *dfdmu0_fp;

  FILE *lfp;
  FILE *efp;

  success = 1;
  dzero            = 0.0;
  num_rxns         = state->number_reactions;
  num_species      = state->nunique_molecules;
  m_r_rt           = state->m_r_rt;
  lfp              = state->lfp;
  ke               = state->ke;
  rke              = state->rke;
  counts           = state->ode_counts;
  conc_to_count    = state->conc_to_count;
  activities       = state->activities;
  molecules_matrix = state->molecules_matrix;
  molecules_ptrs   = molecules_matrix->molecules_ptrs;
  rxn_indices      = molecules_matrix->reaction_indices;
  mcoefficients    = molecules_matrix->coefficients;
  rxn_matrix       = state->reactions_matrix;
  rxn_ptrs         = rxn_matrix->rxn_ptrs;
  molecule_indices = rxn_matrix->molecules_indices;
  rcoefficients    = rxn_matrix->coefficients;
  dfdke_file       = state->dfdke_file;
  dfdmu0_file      = state->dfdmu0_file;
  delta_concs_choice = state->delta_concs_choice;
  dfdke_dfdmu0_work  = state->dfdke_dfdmu0_work;
  /*
    Open the files.
  */
  dfdke_fp = fopen(dfdke_file,"w");
  if (dfdke_fp == NULL) {
    success = 0;
    if (lfp) {
      fprintf(lfp,"compute_dfke_dfdmu0: Error could not open %s\n",
	      dfdke_file);
      fflush(lfp);
    }
  }
  dfdmu0_fp = NULL;
  if (success) {
    state->dfdke_fp = dfdke_fp;
    dfdmu0_fp = fopen(dfdmu0_file,"w");
    if (dfdmu0_fp == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"compute_dfke_dfdmu0: Error could not open %s\n",
		dfdmu0_file);
	fflush(lfp);
      }
    }
  }
  if (success) {
    /*
      Print the dfdke file header.
    */
    fprintf(dfdke_fp,"Partials of concentration derivatives with respect to equilibrium constants\n");
    fprintf(dfdke_fp,"For map of molecule number to name, see .dict file.\n");
    /*
      Print the dfdmu0 file header.
    */
    fprintf(dfdmu0_fp,"Partials of concentration derivatives with respect to chemical potientials\n");
    fprintf(dfdmu0_fp,"For map of molecule number to name, see .dict file.\n");

    
  }
  switch (delta_concs_choice) {
  case 8:
    /*
      Here we have
      
      df_m/dt =  sum  gamma_im (Keq(i) * R'(i)/P_f'''(i) - Keq(i)^(-1)P'(i)/R_r'''(i))
                  i

    */
    telescoping = 0.0;
    forward_piece = dfdke_dfdmu0_work;
    reverse_piece = &forward_piece[num_species];
    mu0_partials  = &reverse_piece[num_species];
    for (i=0;i<num_rxns;i++) {
      pt = 1.0;
      rt = 1.0;
      tr = 1.0;
      tp = 1.0;
      for (j=rxn_ptrs[i];j<rxn_ptrs[i+1];j++) {
        mi = molecule_indices[j];
        /*
        molecule = (struct molecule_struct *)&molecules[mi];
        ci = molecule->c_index;
        compartment = (struct compartment_struct *)&compartments[ci];
        recip_volume       = compartment->recip_volume;
        */
        gamma_im = rcoefficients[j];
        count_mi = counts[mi];
        if (gamma_im < 0.0) {
	  rt = rt * conc_to_pow(count_mi,-gamma_im,telescoping);
	  count_plus = count_mi - gamma_im;
	  tr = tr * conc_to_pow(count_plus,-gamma_im,telescoping);
	  /*
      	  for (k=0;k<(-gamma);k++) {
      	    rt = rt * count_mi;
      	    tr = tr * (count_mi - gamma);
      	  } 
	  */
        } else {
	  if (gamma_im > 0.0) {
	    pt = pt * conc_to_pow(count_mi,gamma_im,telescoping);
	    count_plus = count_mi + gamma_im;
	    tp = tp * conc_to_pow(count_plus,gamma_im,telescoping);
	    /*
	    for (k=0;k<gamma;k++) {
	      pt = pt * count_mi;
	      tp = tp * (count_mi + gamma);
	    }
	    */
	  }
	} 
	forward_piece[i] = (rt/tp);
	reverse_piece[i] = rke[i]* rke[i] * (pt/tr);
	/*
	 */
      } /* end for j */
    }/* end for (i...) */
    /*
      Now loop over the species printing out the partials for each one.
    */
    for (i=0;i<num_species;i++) {
      /*
	zero out a mu0_partials vector.
      */
      vec_set_constant(num_species,mu0_partials,dzero);
      fprintf(dfdke_fp,"molecule %d:  df_%d\n",i,i);
      fprintf(dfdmu0_fp,"molecule %d:  df_%d\n",i,i);
      for (j=molecules_ptrs[i];j<molecules_ptrs[i+1];j++) {
	rj = rxn_indices[j];
	gamma_im = (double)mcoefficients[j];
	dfdke = gamma_im * (forward_piece[rj] + reverse_piece[rj]) * activities[rj];
	if (dfdke != 0) {
	  fprintf(dfdke_fp,"/dke(%d) = %le\n",rj,dfdke);
	  multiplier = dfdke * m_r_rt * ke[rj];
	  for (k=rxn_ptrs[rj];k<rxn_ptrs[rj+1];k++) {
	    mk = molecule_indices[k];
	    gamma_ik = (double)rcoefficients[k];
	    if (gamma_ik != 0) {
	      mu0_partials[mk] += multiplier * gamma_ik;
	    }
	  } /* end for (k...) */
	}
      } /* end for (j...) */
      for (k=0;k<num_species;k++) {
	if (mu0_partials[k] != dzero) {
	  fprintf(dfdmu0_fp,"/dmu(%d) = %le\n",k,mu0_partials[k]);
	}
      } /* end for (k...) */
    } /* end for (i...) */
    break;
  } /* end switch(delta_concs_choice) */
  if (dfdke_fp) {
    fclose(dfdke_fp);
  }
  if (dfdmu0_fp) {
    fclose(dfdmu0_fp);
  }
  return(success);
}
