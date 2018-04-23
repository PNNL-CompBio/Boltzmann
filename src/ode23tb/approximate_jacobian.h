#ifndef _APPROXIMATE_JACOBIAN_H_
#define _APPROXIMATE_JACOBIAN_H_ 1
extern int approximate_jacobian(struct state_struct *state, 
				double *concs,
				double *delta_concs,
				double t,
				int choice);
#endif
