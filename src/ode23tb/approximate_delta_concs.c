#include "boltzmann_structs.h"

#include "lr_approximate_delta_concs.h"
#include "lr1_approximate_delta_concs.h"
#include "lr2_approximate_delta_concs.h"
#include "lr3_approximate_delta_concs.h"
#include "lr4_approximate_delta_concs.h"
#include "lr5_approximate_delta_concs.h"
#include "lr6_approximate_delta_concs.h"
#include "lr7_approximate_delta_concs.h"
#include "lr8_approximate_delta_concs.h"
#include "lr9_approximate_delta_concs.h"
#include "lr10_approximate_delta_concs.h"
#include "lr11_approximate_delta_concs.h"
#include "lr12_approximate_delta_concs.h"
#include "lr13_approximate_delta_concs.h"
#include "lr14_approximate_delta_concs.h"

#include "approximate_delta_concs.h"

int approximate_delta_concs(struct state_struct *state, double *concs,
			    double *flux, int choice) {
  /*
    Compute approximations to concentartion changes wrt time
    0 for lr_approximate_delta_concs, based on likelihood ratios.
    1 for lr1_approximate_delta_concs, based on likelihood ratios.
    2 for lr2_approximate_delta_concs, based on likelihood ratios.
    3 for lr3_approximate_delta_concs, based on likelihood ratios.


    Called by: ode23tb, num_jac_col, ode_it_solve
    Calls      lr_approximate_delta_concs,
	       

                                TMF
    state                       *SI   Boltzmant state structure.
                                      uses number_reactions,
				           unique_moleules,
                                           molecules_matrix,
					   deriv_thresh,
					   and lfp,
				      

    concs			D1I   molecule concentrations vector of length 
                                      nunique_moleucles

    flux                        D1O   vector of length  unique_molecules
                                      of concentration change per unit time.
				      Set by this routine.

    choice                      I0I   0 for lr_approximate_delta_concs


  */
  double deriv_thresh;
  int success;
  int i;
  int ny;
  int padi;
  FILE *lfp;
  FILE *efp;
  lfp = state->lfp;
  success = 1;
  switch(choice) {
  case 0:
    /*
      Default Likelihood ratio approximation.
    */
    success = lr_approximate_delta_concs(state,
					 concs,
					 flux,
					 choice);
    break;
  case 1:
    success = lr1_approximate_delta_concs(state,
					  concs,
					  flux,
					  choice);
    break;
  case 2:
    success = lr2_approximate_delta_concs(state,
					  concs,
					  flux,
					  choice);
    break;
  case 3:
    success = lr3_approximate_delta_concs(state,
					  concs,
					  flux,
					  choice);
    break;
  case 4:
    success = lr4_approximate_delta_concs(state,
					  concs,
					  flux,
					  choice);

    break;
  case 5:
    success = lr5_approximate_delta_concs(state,
					  concs,
					  flux, 
					  choice);

    break;
  case 6:
    success = lr6_approximate_delta_concs(state,
					  concs,
					  flux, 
					  choice);

    break;
  case 7:
    success = lr7_approximate_delta_concs(state,
					  concs,
					  flux, 
					  choice);

    break;
  case 8:
    success = lr8_approximate_delta_concs(state,
					  concs,
					  flux, 
					  choice);

    break;
  case 9:
    success = lr9_approximate_delta_concs(state,
					  concs,
					  flux, 
					  choice);

    break;
  case 10:
    success = lr10_approximate_delta_concs(state,
					   concs,
					   flux, 
					   choice);
    break;
  case 11:
    success = lr11_approximate_delta_concs(state,
					   concs,
					   flux, 
					   choice);
    break;
  case 12:
    success = lr12_approximate_delta_concs(state,
					   concs,
					   flux, 
					   choice);
    break;
  case 13:
    success = lr13_approximate_delta_concs(state,
					   concs,
					   flux, 
					   choice);
    break;
  case 14:
    success = lr14_approximate_delta_concs(state,
					   concs,
					   flux, 
					   choice);
    break;
    /*
      Use debugging flavor with kinetic rate constants for coupledenzyme.in
    */
  default :
    success = 0;
    if (lfp) {
      fprintf(lfp,"approximate_delta_concs: invalid delta_concs_choice was %d\n",choice);
      fflush(lfp);
    }
  } /* end switch (choice) */
  /*
    Post process the flux vector to zero out fluxes less than deriv_thresh;
  */
  ny = state->nunique_molecules;
  deriv_thresh = state->deriv_thresh;
  for (i=0;i<ny;i++) {
    if (fabs(flux[i]) < deriv_thresh) {
      flux[i] = 0.0;
    }
  }
  return(success);
}