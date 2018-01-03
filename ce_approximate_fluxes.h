#ifndef _CE_APPROXIMATE_FLUXES_H_
#define _CE_APPROXIMATE_FLUXES_H_ 1
extern int ce_approximate_fluxes(struct state_struct *state, double *counts,
				 double *forward_rxn_likelihoods,
				 double *reverse_rxn_likelihoods, 
				 double *flux, double multiplier,
				 int base_rxn);
#endif
