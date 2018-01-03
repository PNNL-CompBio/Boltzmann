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
  struct rxn_matrix_struct *reactions_matrix; 
  int64_t *molecules_indices;
  int64_t *coefficients;
  int64_t *rxn_ptrs;
  double *kf_rel;
  double *kr_rel;
  double *delta_g0_tfs;
  double forward_g0_sum;
  double reverse_g0_sum;
  double recip_rt;
  int64_t j;
  int64_t coefficientj;
  int64_t moleculej;
  int nrxns;
  int success;
  int i;
  int padi;

  nrxns             = (int)state->number_reactions;
  delta_g0_tfs      = state->molecule_dg0tfs;
  reactions_matrix  = state->reactions_matrix;
  molecules_indices = reactions_matrix->molecules_indices;
  coefficients      = reactions_matrix->coefficients;
  rxn_ptrs          = reactions_matrix->rxn_ptrs;
  recip_rt          = 1.0/state->rt;
  kf_rel            = state->kf_rel;
  kr_rel            = state->kr_rel;
  success = 1;
  for (i=0;i < nrxns;i++) {
    forward_g0_sum = 0.0;
    reverse_g0_sum = 0.0;
    for (j=rxn_ptrs[i];j<rxn_ptrs[i+1];j++) {
      moleculej = molecules_indices[j];
      coefficientj = coefficients[j];
      if (coefficientj < 0) {
	forward_g0_sum -= (coefficientj * delta_g0_tfs[moleculej]);
      } else {
	if (coefficientj > 0) {
	  reverse_g0_sum += (coefficientj * delta_g0_tfs[moleculej]);
	}
      }
    }
    kf_rel[i] = exp(forward_g0_sum * recip_rt);
    kr_rel[i] = exp(reverse_g0_sum * recip_rt);
  }
  return(success);
}
  
