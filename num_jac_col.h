#ifndef _NUM_JAC_COL_H_
#define _NUM_JAC_COL_H_ 1
extern int num_jac_col(struct state_struct *state,
	   	       int ny, int j,
		       int *rowmax_p,
		       double *y,
		       double *f,
		       double *delj_p,
		       double threshj,
		       double *fdel,
		       double *fdiff,
		       double *dfdy_colj,
		       double *absfvaluerm_p,
		       double *absfdelrm_p,
		       double *absfdiffmax_p,
		       double *infnormdfdy_colj_p);
#endif
