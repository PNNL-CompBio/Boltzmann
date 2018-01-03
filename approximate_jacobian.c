#include "boltzmann_structs.h"
#include "lr8_approximate_jacobian.h"

int approximate_jacobian(struct state_struct *state, 
			 double *concs,
			 double *delta_concs,
			 double t,
			 int choice) {
  /*
    Compute approximations to the jacobian of theconcentration changes. wrt time,   Based on thermodynamics formulation for concentraion rate
    changes     using counts to compute tr, tp, pt, rt instead of concs.
    This routine uses the following state fields,
          rxns_matrix, 
	  unique_molecules, 
	  sorted_molecules,
	  ke, 
	  rke,
	  cvodes_params,
     It sets the following fields of cvodes_params,
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
  */
  int success;
  int padi;
  success = 1;
  switch (choice) {
  case 0:
    /*
      This should be a numerical jacobian approximation.
    */
    success = 0;
    break;
  case 8:
    success = lr8_approximate_jacobian(state,concs,delta_concs,t,choice);
    break;
  }
  return(success);
}
