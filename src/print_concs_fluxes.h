#ifndef _PRINT_CONCS_FLUXES_H_
#define _PRINT_CONCS_FLUXES_H_ 1
extern int print_concs_fluxes(struct state_struct *state,int ny,
			      double *fluxes, 
			      double *concs, 
			      double time, double h,
			      int step, int origin);
#endif
