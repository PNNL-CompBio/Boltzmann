#ifndef _ODE_NUM_JAC_H_
#define _ODE_NUM_JAC_H_ 1
extern int ode_num_jac(struct state_struct *state,
		       int    first_time,
		       double *dfdy, 
		       double t, 
		       double *y, 
		       double *f, 
		       double *fac,
		       double *thresh,
		       double *fdel, 
		       double *fdiff, 
		       double *dfdy_tmp, 
 		       int64_t *nf);
#endif
