#ifndef _ODE23TB_INIT_H_H_
#define _ODE23TB_INIT_H_H_ 1
extern void ode23tb_init_h(int normcontrol, 
			   int ny, 
			   double normy, 
			   double normyp,
			   double recip_70p,
			   double recip_cube_root_rtol, 
			   double htspan, 
			   double hmin, 
			   double hmax, 
			   double tdir,
			   double *y, 
			   double *yp,
			   double *wt,
			   double *h_p,
			   double *absh_p);
#endif
