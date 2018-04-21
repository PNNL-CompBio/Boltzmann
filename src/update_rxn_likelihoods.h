#ifndef _UPDATE_RXN_LIKELIHOODS_H_
#define _UPDATE_RXN_LIKELIHOODS_H_ 1
extern int update_rxn_likelihoods(struct state_struct *state, double *counts,
				  double *forward_rxn_likelihood,
				  double *reverse_rxn_likelihood);
#endif
