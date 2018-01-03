#include "boltzmann_structs.h"
#include "ode_num_jac.h"
#include "lr8_approximate_jacobian.h"
#include "boltzmann_dense_to_sparse.h"
#include "boltzmann_sparse_to_dense.h"
#include "approximate_jacobian.h"
int approximate_jacobian(struct state_struct *state, 
			 double *concs,
			 double *delta_concs,
			 double t,
			 int choice) {
  /*
    Compute approximations to the jacobian of the concentration changes. 
    wrt time,   Based on thermodynamics formulation for concentraion rate
    changes     using counts to compute tr, tp, pt, rt instead of concs.
    This routine uses the following state fields,
          rxns_matrix, 
	  unique_molecules, 
	  sorted_molecules,
	  ode_num_jac_first_call,
	  ke, 
	  rke,
	  ode23tb_params,
	  cvodes_params,
	  
     It sets the following fields of cvodes_params,
     if choice is not 0 or solver choice is not 0.
	  prec_row (= rfc),
	  drfc,
	  dfdy_a,
	  dfdy_ia,
	  dfdy_ja,
	  dfdy_at,
	  dfdy_iat,
	  dfdy_jat,
	  column_index,
	  column_mask

     It sets the following fields of ode23tb_params,
     if choice is 0 or solver choice is 0.
          dfdy,
	  fac,
	  thresh,
	  fdel,
	  fiff,
	  dfdy_tmp,
	  


     We might want it to have persistent scratch spaces dfdy(ny,ny),
     fac(ny), thresh(ny), fdel(ny), fdiff(ny), dfdy_tmp(ny) and 
     scalar ode_num_jac_first_time, these used for call to ode_num_jac.


     Called by: boltzmann_cvodes_jtimes
     Calls:
  */
  struct ode23tb_params_struct *ode23tb_params;
  struct cvodes_params_struct  *cvodes_params;
  double *dfdy;
  double *fac;
  double *thresh;
  double *fdel;
  double *fdiff;
  double *dfdy_tmp;
  double *dfdy_a;
  int  *dfdy_ia;
  int  *dfdy_ja;
  int64_t nf;
  int success;
  int first_time;
  int ode_solver_choice;
  int ny;
  success = 1;
  ode23tb_params = state->ode23tb_params;
  cvodes_params  = state->cvodes_params;
  ode_solver_choice = state->ode_solver_choice;
  ny                = (int)state->nunique_molecules;
  switch (choice) {
  case 0:
    /*
      This should be a numerical jacobian approximation.
    */
    nf       = ode23tb_params->nf;
    dfdy     = ode23tb_params->dfdy;
    fac      = ode23tb_params->fac;
    thresh   = ode23tb_params->thresh;
    fdel     = ode23tb_params->fdel;
    fdiff    = ode23tb_params->fdiff;
    dfdy_tmp = ode23tb_params->dfdy_tmp;
    first_time = ode23tb_params->num_jac_first_time;
    ode_num_jac(state,first_time,
		dfdy,t,concs,delta_concs,
		fac,thresh,fdel,fdiff,dfdy_tmp,&nf);
    ode23tb_params->nf = nf;
    ode23tb_params->num_jac_first_time = 0;
    /*
      Now if the ode choice is not ode23tb we need to 
      convert dfdy to dfdy_a, dfdy_ia, dfdy_ja for 
      use by cvodes.
    */
    if (ode_solver_choice != 0) {
      dfdy_a = cvodes_params->dfdy_a;
      dfdy_ia = cvodes_params->dfdy_ia;
      dfdy_ja = cvodes_params->dfdy_ja;
      boltzmann_dense_to_sparse(ny,dfdy,dfdy_a,dfdy_ia,dfdy_ja);
    }
    break;
  case 8:
    success = lr8_approximate_jacobian(state,concs,delta_concs,t,choice);
    /*
      Now if the ode choice is not cvodes we need to convert
      dfdy_a, dfdy_ia, dfdy_ja to dfdy for ode23tb.
    */
    if (ode_solver_choice == 0) {
      dfdy_a  = cvodes_params->dfdy_a;
      dfdy_ia = cvodes_params->dfdy_ia;
      dfdy_ja = cvodes_params->dfdy_ja;
      dfdy    = ode23tb_params->dfdy;
      boltzmann_sparse_to_dense(ny,dfdy_a,dfdy_ia,dfdy_ja,dfdy);
    }
    break;
  } /* end switch (choicee) */
  return(success);
}
