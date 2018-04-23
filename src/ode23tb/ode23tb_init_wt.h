#ifndef _ODE23TB_INIT_WT_H_
#define _ODE23TB_INIT_WT_H_ 1
extern void ode23tb_init_wt(int normcontrol, 
			    int ny, 
			    double normy, 
			    double threshold,
			    double *y, 
			    double *wt);
#endif
