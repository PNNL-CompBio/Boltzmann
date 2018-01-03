#include "boltzmann_structs.h"
#include "init_relative_rates.h"
int init_relative_rates(struct state_struct *state) {
  /*
    Compute the relative reaction rates, forward, kf_rel, and reverse
    kr_rel, from deltag0's and reaction matrix, for use in 
    lr6_approximate_delta_concs.
    This is "computational alchemy" - ask Bill Cannon for reference.
    kf_rel and kr_rel are allocated in alloc7 routine.

    Called by: deq_run.
    Calls:
  */
  struct reactions_matrix_struct *reactions_matrix; 
  int64_t *molecules_indices;
  int64_t *coefficients;
  int64_t *rxn_ptrs;
  double *log_kf_rel;
  double *log_kr_rel;
  double *delta_g0_tfs;
  double forward_g0_sum;
  double reverse_g0_sum;
  double m_recip_rt;
  double scaled_fg0s;
  double scaled_rg0s;
  /*
  double g0_sum_min;
  double g0_sum_max;
  double g0_range_mid;
  double adj_f_g0;
  double adj_r_g0;
  */
  int64_t j;
  int64_t coefficientj;
  int64_t moleculej;
  int nrxns;
  int success;
  int i;
  int padi;
  FILE *lfp;
  FILE *efp;

  nrxns             = (int)state->number_reactions;
  delta_g0_tfs      = state->molecule_dg0tfs;
  reactions_matrix  = state->reactions_matrix;
  molecules_indices = reactions_matrix->molecules_indices;
  coefficients      = reactions_matrix->coefficients;
  rxn_ptrs          = reactions_matrix->rxn_ptrs;
  m_recip_rt        = -1.0/state->rt;
  log_kf_rel        = state->log_kf_rel;
  log_kr_rel        = state->log_kr_rel;
  lfp               = state->lfp;
  success = 1;
  for (i=0;i < nrxns;i++) {
    forward_g0_sum = 0.0;
    reverse_g0_sum = 0.0;
    for (j=rxn_ptrs[i];j<rxn_ptrs[i+1];j++) {
      moleculej = molecules_indices[j];
      coefficientj = coefficients[j];
      if (coefficientj < 0) {
	forward_g0_sum += (coefficientj * delta_g0_tfs[moleculej]);
      } else {
	if (coefficientj > 0) {
	  reverse_g0_sum -= (coefficientj * delta_g0_tfs[moleculej]);
	}
      }
    }
    /*
    kf_rel[i] = exp(forward_g0_sum * m_recip_rt);
    kr_rel[i] = exp(reverse_g0_sum * m_recip_rt);
    */
    scaled_fg0s = forward_g0_sum * m_recip_rt;
    scaled_rg0s = reverse_g0_sum * m_recip_rt;
    log_kf_rel[i] = scaled_fg0s;
    log_kr_rel[i] = scaled_rg0s;
  }
  /*
    Scale the g0's to be centered about zero in the union of their ranges.
    could look like
    g0_sum_min = kf_rel[0];
    min_cand   = kr_rel[0];
    max_cand   = min_cand + g0_sum_min;
    new_g0_min = min_cand < g0_sum_min;
    new_mask   = 0 - ((int64_t)new_g0_min);
    cand_delta = (min_cand - g0_sum_min); 
    // probably need to do this op with a union no way to type cast.
    memcopy(maskable_delta,cand_delta,8);
    maskable_delta  =  maskable_delta & new_mask;
    memcopy(cand_delta,maskable_delta,8);
    g0_sum_min = g0_sum_min + cand_delta;
    g0_sum_max = max_cand - g0_sum_min;
    for (i=1;i<nrxns;i++) {
       min_cand = kf_rel[i];
       max_cand = kr_rel[i];
       sum_cand = max_cand + min_cand;
       swap_cand = max_cand < min_cand;
       swap_mask = 0 - ((int64_t)swap_cand;
       cand_delta = max_cand - min_cand;
       memcpy(maskable_delta,delta_cand,8);
       maskable_delta = maskable_delta & swap_mask;
       memcpy(maskable_delta,delta_cand,8);
       min_cand = min_cand + delta_cand;
       max_cand = max_cand - delta_cand;

       new_g0_min = min_cand < g0_sum_min;
       new_mask   = 0 - ((int64_t)new_g0_min);
       cand_delta = min_cand - g0_sum_min;
       memcpy(maskable_delta,cand_delta,8);
       maskable_delta  =  maskable_delta & new_mask;
       memcpy(cand_delta,maskable_delta,8);
       g0_sum_min = g0_sum_min + cand_delta;
       new_g0_max = max_cand > g0_sum_max;
       new_mask = 0 - ((int64_t)new_g0_max);
       cand_delta = max_cand - g0_sum_max;
       memcpy(maskable_delta,cand_delta,8);
       maskable_delta  =  maskable_delta & new_mask;
       memcpy(cand_delta,maskable_delta,8);
       g0_sum_max = g0_sum_max + cand_delta;

       Nice idea but it doesn't work as 
       the range of kf/kr pairs is still to broad.
  if (kf_rel[0] < kr_rel[0]) {
    g0_sum_min = kf_rel[0];
    g0_sum_max = kr_rel[0];
  } else {
    g0_sum_min = kr_rel[0];
    g0_sum_max = kf_rel[0];
  }
  for (i=1;i<nrxns;i++) {
    if (kf_rel[i] < g0_sum_min) {
      g0_sum_min = kf_rel[i];
    } else {
      if (kf_rel[i] > g0_sum_max) {
	g0_sum_max = kf_rel[i];
      }
    }
    if (kr_rel[i] < g0_sum_min) {
      g0_sum_min = kr_rel[i];
    } else {
      if (kr_rel[i] > g0_sum_max) {
	g0_sum_max = kr_rel[i];
      }
    }
  }
  g0_range_mid = 0.5 * (g0_sum_max + g0_sum_min);
  for (i=0;i<nrxns;i++) {
    adj_f_g0 = kf_rel[i] - g0_range_mid;
    kf_rel[i] = exp(adj_f_g0 * m_recip_rt);
    adj_r_g0 = kr_rel[i] - g0_range_mid;
    kr_rel[i] = exp(adj_f_g0 * m_recip_rt);
  }

  */

#ifdef DBG_INIT_RELATIVE_RATES
  lfp = state->lfp;
  if (lfp) {
    for (i=0;i<nrxns;i++) {
      fprintf(lfp,"rxn %d, log_kf_rel = %le, log_kr_rel = %le\n",
	      i,log_kf_rel[i],log_kr_rel[i]);
    }
    fflush(lfp);
  }
#endif
  return(success);
}
  
