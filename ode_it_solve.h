#ifndef _ODE_IT_SOLVE_H_
#define _ODE_IT_SOLVE_H_ 1
extern int ode_it_solve(struct state_struct *state,
			double t,
			double *y,
			double *z,
			double *del,
			double *znew,
			double *y_count,
			double *forward_rxn_likelihoods,
			double *reverse_rxn_likelihoods,
			double d,
			double h,
			double rtol,
			double wt,
			double *rate,
			int *iter_count);
#endif
