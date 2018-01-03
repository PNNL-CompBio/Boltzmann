#ifndef _BOLTZMANN_CVODES_PSOLVE_H_
#define _BOLTZMANN_CVODES_PSOLVE_H_ 1
extern int boltzmann_cvodes_psolve(double t,
				   N_Vector y,
				   N_Vector fy,
				   N_Vector r,
				   N_Vector z,
				   double gamma,
				   double delta,
				   int lr,
				   void *user_data,
				   N_Vector tmp);
#endif
