#ifndef _ODE_SOLVER_H_
#define _ODE_SOLVER_H_ 1
extern int ode_solver (struct state_struct *state, double *counts,
		       double htry, int nonnegative, int normcontrol,
		       int print_concs, int choice);
#endif
