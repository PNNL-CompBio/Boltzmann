#ifndef _PRINT_CONCS_GRAD_H_
#define _PRINT_CONCS_GRAD_H_ 1
extern int print_concs_grad(struct state_struct *state,int ny,
			      double *grad, 
			      double *concs, 
			      double time, double h,
			      int step, int origin);
#endif
