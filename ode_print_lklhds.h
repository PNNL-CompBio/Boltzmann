#ifndef _ODE_PRINT_LKLHDS_H_
#define _ODE_PRINT_LKLHDS_H_ 1
extern void ode_print_lklhds(struct state_struct *state,
			     double t,
			     double *forward_rxn_likelihoods,
			     double *reverse_rxn_likelihoods);
#endif
