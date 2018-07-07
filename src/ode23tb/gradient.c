#include "boltzmann_structs.h"

#include "lr0_gradient.h"
#include "lr1_gradient.h"
#include "lr2_gradient.h"
#include "lr3_gradient.h"
#include "lr4_gradient.h"
#include "lr5_gradient.h"
#include "lr6_gradient.h"
#include "lr7_gradient.h"
#include "lr8_gradient.h"
#include "lr9_gradient.h"
#include "lr10_gradient.h"
#include "lr11_gradient.h"
#include "lr12_gradient.h"
#include "lr13_gradient.h"
#include "lr14_gradient.h"

#include "gradient.h"

int gradient(struct state_struct *state, double *concs,
			    double *flux, int choice) {
  /*
    Compute approximations to concentartion changes wrt time
    0 for lr_gradient, based on likelihood ratios.
    1 for lr1_gradient, based on likelihood ratios.
    2 for lr2_gradient, based on likelihood ratios.
    3 for lr3_gradient, based on likelihood ratios.


    Called by: ode23tb, num_jac_col, ode_it_solve
    Calls      lr_gradient,
	       

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

    choice                      I0I   0 for lr_gradient


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
    success = lr0_gradient(state,concs,flux,choice);
    break;
  case 1:
    success = lr1_gradient(state,concs,flux,choice);
    break;
  case 2:
    success = lr2_gradient(state,concs,flux,choice);
    break;
  case 3:
    success = lr3_gradient(state,concs,flux,choice);
    break;
  case 4:
    success = lr4_gradient(state,concs,flux,choice);
    break;
  case 5:
    success = lr5_gradient(state,concs,flux,choice);
    break;
  case 6:
    success = lr6_gradient(state,concs,flux,choice);
    break;
  case 7:
    success = lr7_gradient(state,concs,flux,choice);
    break;
  case 8:
    success = lr8_gradient(state,concs,flux,choice);
    break;
  case 9:
    success = lr9_gradient(state,concs,flux,choice);
    break;
  case 10:
    success = lr10_gradient(state,concs,flux,choice);
    break;
  case 11:
    success = lr11_gradient(state,concs,flux,choice);
    break;
  case 12:
    success = lr12_gradient(state,concs,flux,choice);
    break;
  case 13:
    success = lr13_gradient(state,concs,flux,choice);
    break;
  case 14:
    success = lr14_gradient(state,concs,flux,choice);
    break;
    /*
      Use debugging flavor with kinetic rate constants for coupledenzyme.in
    */
  default :
    success = 0;
    if (lfp) {
      fprintf(lfp,"gradient: invalid gradient_choice was %d\n",choice);
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
