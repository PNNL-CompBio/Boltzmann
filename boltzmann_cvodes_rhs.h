#ifndef _BOLTZMANN_CVODES_RHS_H_
#define _BOLTZMANN_CVODES_RHS_H_ 1
extern int boltzmann_cvodes_rhs(double t, N_Vector y, N_Vector y_dot,
				void *user_data);
#endif
