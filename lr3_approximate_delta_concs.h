#ifndef _LR3_APPROXIMATE_DELTA_CONCS_H_
#define _LR3_APPROXIMATE_DELTA_CONCS_H_ 1
extern int lr3_approximate_delta_concs(struct state_struct *state, 
				       double *counts,
				       double *forward_rxn_likelihoods,
				       double *reverse_rxn_likelihoods, 
				       double *flux, double flux_scaling,
				       int base_rxn, int choice);
#endif
