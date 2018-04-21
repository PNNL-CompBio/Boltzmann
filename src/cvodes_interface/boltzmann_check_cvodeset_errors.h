#ifndef _BOLTZMANN_CHECK_CVODESET_ERRORS_H_
#define _BOLTZMANN_CHECK_CVODESET_ERRORS_H_ 1
extern int boltzmann_check_cvodeset_errors(int flag,
					   void *cvode_mem,
					   struct state_struct *state,
					   const char *routine);
#endif
