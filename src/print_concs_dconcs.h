#ifndef _PRINT_CONCS_DCONCS_H_
#define _PRINT_CONCS_DCONCS_H_ 1
extern int print_concs_dconcs(struct state_struct *state,int ny,
			      double *dconcs, 
			      double *concs, 
			      double time, double h,
			      int step, int origin);
#endif
