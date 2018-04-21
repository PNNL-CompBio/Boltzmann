#ifndef _ODE23TB_PARAMS_STRUCT_H_
#define _ODE23TB_PARAMS_STRUCT_H_ 1
struct ode23tb_params_struct {
  /*
    Workspace for jacobian computation, with numerical approximation.
    dfdy, the jacobian dense matrix, length is ny * ny.
  */
  double *dfdy; 
  /*
    History mechanism, length ny.
  */
  double *fac;
  /*
    Threshold vector. length ny
  */
  double *thresh;
  /*
    fdel, length ny
  */
  double *fdel;
  /*
    fdiff, length ny
  */
  double *fdiff;
  /*
    dfdy_tmp, length ny
  */
  double *dfdy_tmp;
  /*
    Initial step size.
  */
  double htry;
  /*
    threshold, default value 1.e-3
  */
  double threshold; 
  /*
    Numerical Jacobian threshold; default value 1.e-6
  */
  double njthreshold;
  /*
    Non-negative flag 1 = require non-negative result, 0 = no requirement
  */
  int nonnegative;
  /*
    Norm control. 0 for infinity, 1 for 2-norm, 
  */
  int normcontrol;
  /*
    First time flag for ode_num_jac. 1 for yes, 0 for no.
  */
  int num_jac_first_time;
  /*
    number of function calls.
  */
  int nf;
}
;
#endif
