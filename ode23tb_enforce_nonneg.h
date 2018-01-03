#ifndef _ODE23TB_ENFORCE_NONNEG_H_
#define _ODE23TB_ENFORCE_NONNEG_H_ 1
extern void ode23tb_enforce_nonneg(int ny, int normcontrol, double *ynew,
				   double *znew, double *normynew_p);
#endif
