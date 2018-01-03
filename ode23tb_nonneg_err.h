#ifndef _ODE23TB_NONNEG_ERR_H_
#define _ODE23TB_NONNEG_ERR_H_ 1
extern void ode23tb_nonneg_err(int ny, 
			       double *ynew,
			       double rtol, 
			       double ernn_scale,
			       double *err_p,
			       int *nnrejectstep_p);
#endif
