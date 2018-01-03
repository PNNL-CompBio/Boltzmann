#ifndef _BOLTZMANN_CVODES_PSETUP_H_
#define _BOLTZMANN_CVODES_PSETUP_H_ 1
extern int boltzmann_cvodes_psetup(double t, 
				   N_Vector y,
				   N_Vector fy,
				   int jok,
				   int *jcurptr,
				   double gamma,
				   void *user_data,
				   N_Vector tmp1,
				   N_Vector tmp2,
				   N_Vector tmp3);
#endif
