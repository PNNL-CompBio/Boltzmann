#ifndef _ODE_IT_SOLVE_H_
#define _ODE_IT_SOLVE_H_ 1
extern int ode_it_solve(struct state_struct *state,
			double *miter,
			int    *ipivot,
			double t,
			double *y,
			double *z,
			double *del,
			double *rhs,
			double *norm_space,
			double d,
			double h,
			double rtol,
			double *wt,
			double *rate,
			int *iter_count_p);
#endif
