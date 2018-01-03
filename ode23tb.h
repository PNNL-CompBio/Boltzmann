#ifndef _ODE23TB_H_
#define _ODE23TB_H_ 1
extern int ode23tb (struct state_struct *state, double *counts,
		    double htry, int nonnegative,
		    int normcontrol);
#endif
