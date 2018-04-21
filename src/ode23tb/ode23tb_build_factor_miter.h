#ifndef _ODE_BUILD_FACTOR_MITER_H_
#define _ODE_BUILD_FACTOR_MITER_H_ 1
extern int ode23tb_build_factor_miter(int ny,
				      int nysq,
				      double d,
				      double h,
				      double *dfdy,
				      double *miter,
				      int    *ipivot,
				      int    *info_p,
				      FILE   *lfp);
#endif
