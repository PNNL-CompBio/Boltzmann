#ifndef _BOLTZMANN_CVODES_JTIMES_H_
#define _BOLTZMANN_CVODES_JTIMES_H_ 1
extern int boltzmann_cvodes_jtimes(N_Vector v,
				   N_Vector jv,
				   double t,
				   N_Vector y,
				   N_Vector fy,
				   void *user_data,
				   N_Vector tmp);
#endif
