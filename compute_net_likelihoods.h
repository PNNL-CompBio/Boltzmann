#ifndef _COMPUTE_NET_LIKELIHOODS_H_
#define _COMPUTE_NET_LIKELIHOODS_H_ 1
extern int compute_net_likelihoods(struct state_struct *state,
				   double *forward_rxn_likeklihood,
				   double *reverse_rxn_likeklihood,
				   double *net_likelihood);
#endif
