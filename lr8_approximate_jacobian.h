#ifndef _LR8_APPROXIMATE_JACOBIAN_H_
#define _LR8_APPROXIMATE_JACOBIAN_H_
extern int lr8_approximate_jacobian(struct state_struct *state, 
				    double *concs,
				    double *delta_concs,
				    double t,
				    int choice);
#endif
