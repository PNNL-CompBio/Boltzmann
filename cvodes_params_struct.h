#ifndef _CVODES_PARAMS_STRUCT_H_
#define _CVODES_PARAMS_STRUCT_H_ 1
struct cvodes_params_struct {
  /*
    Workspace created by CVodeCreate call.
  */
  void *cvode_mem;

  /*
    Pointer to workspace for jacobian setup, factorization and application
    routines.
  */
  double *drfc;
  double *dfdy_a;
  double *dfdy_at;
  double *miter_m;
  double *prec_l;
  double *prec_u;
  double *prec_row;
  double *recip_diag_u;
  double *frow;
  double *srow;
  int    *dfdy_ia;
  int    *dfdy_ja;
  int    *dfdy_iat;
  int    *dfdy_jat;
  int    *miter_im;
  int    *miter_jm;
  int    *prec_il;
  int    *prec_jl;
  int    *prec_iu;
  int    *prec_ju;
  int    *lindex;
  int    *uindex;
  int    *column_mask;
  int    *sindex;
  /*
    Relative tolerance.
  */
  double reltol;
  /*
    Absolute tolerance.
  */
  double abstol;
  /*
    Initial step size.
  */
  double hin;

  /*
    Minium step size.
  */
  double hmin;
  
  /*
    Maximum step size.
  */
  double hmax;

  /*
    Stop time.
  */
  double tstop;
  /*
    Coefficient in the nonlinear convergence trest.
  */
  double nlscoef;
  /*
    Ratio between linear and nonlinear tolerances, default .05
  */
  double eplifac;
  /*
    The linear multistep indicator, one of CV_ADAMS, or CV_BDF.
    As the boltzmann equations are likely to be stiff use CV_BDF
  */
  int linear_multistep_method; 
  /*
    The iterative method indicator, one of CV_NEWTON or CV_FUNCTIONAL,
    Recommended choice is CV_NEWTON.
  */
  int iterative_method;

  /*
    Linear solver choice
  */
  int linear_solver_method;
  /*
    Max number of direction (Krylov) vectors for GMRES, 
    default for cvodes is 5.
    I recommend 30 djb.
  */
  int maxl;

  /*
    Upper and lower halfbandwidths for the band solver.
  */
  int mupper;
  int mlower;
  

  /*
    nnz, number of nonzero entries in Jacobian
  */
  int nnz;
  int sparsetype; /* one of CSC_MAT or CSR_MAT */


  /*
    Preconditioner type. One of PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH
    Recommend PREC_LEFT
  */
  int pretype;
  /*
    Number of steps you want cvode to output data for with itask = CVNORMAL
    -1 for using itask = CV_ONE_STEP.
  */
  int num_cvode_steps;

  /*
    Maximum number of internal steps, default value to be 500.
  */
  int mxsteps;
  /*
    Maximum number of error test failures per internal CVode step.
    Default value is 7.
  */
  int mxnef;
  
  /*
    Maximun number of convergent test failures.
  */
  int mxncf;
  /*
    Maximum order of the linear multistep method.
  */
  int max_ord;

  /*
    Maximum order for Adams_Moulton method = 12
  */
  int adams_q_max;
  /*
    Maximum order for the Backward Differention (BDF) method, = 5;
  */
  int bdf_q_max;

  /*
    Maximum number of warnings for zero step size h. Default 10.
  */
  int mxhnil;
  /*
    Use stablity detection limit. 1 for yes, 0 for no.
  */
  int use_stab_lim_det;

  /*
    Maximum number of error test failures. Default 7.
  */
  int maxnef;
  /*
    Maximum number of nonlinear iteratons. Default 3.
  */
  int maxcor;

  /*
    Maximum number of convergence failures default 10.
  */
  int maxncf;
  /*
    Root direction.
  */
  int root_direction;
  
  /*
    Disable rootfinding warnings flag.
  */
  int no_inactive_root_warn;
  /*
    sparse matrix ordering alogroithm.
    not used.
  */
  int sparse_matrix_ordering_alg;
  /*
    Gram-Schmidt orthogonalization, default modified, MODIFIED_GS
  */
  int gstype;
  /*
    cvodes rhs choice (function selector).
  */
  int cvodes_rhs_choice;

  /*
    cvodes jacobian multplier choice.
  */
  int cvodes_jtimes_choice;
  /*
    cvodes preconditioner choice.
  */
  int cvodes_prec_choice;

  /*
    Maximum length of miter_m, and miter_jm
  */
  int nnzm;
  /*
    Maximum lenght of prec_l and prec_jl
  */
  int nnzl;

  /*
    Maximum lenght of prec_u and prec_ju
  */
  int nnzu;
  /*
    Preconditioner fill level.
  */
  int prec_fill;
}
;
#endif
