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
#include "ce_approximate_delta_concs.h"

#include "approximate_delta_concs.h"

int approximate_delta_concs(struct state_struct *state, double *concs,
			    double *flux, int choice) {
  /*
    Compute approximations to concentartion changes wrt time
    0 for lr_approximate_delta_concs, based on likelihood rations.
    1 for lr1_approximate_delta_concs, based on likelihood rations.
    2 for lr2_approximate_delta_concs, based on likelihood rations.
    3 for lr3_approximate_delta_concs, based on likelihood rations.

    42 for ce_approximate_delta_concs, for debugging only with
    the reaction rates for the coupledenzyme.in file

    Called by: ode23tb, num_jac_col, ode_it_solve
    Calls      lr_approximate_delta_concs,
               ce_approximate_delta_concs,
	       

                                TMF
    state                       *SI   Boltzmant state structure.
                                      uses number_reactions,
				           unique_moleules,
                                           molecules_matrix,
					   and lfp,
				      

    concs			D1I   molecule concentrations vector of length 
                                      nunique_moleucles

    flux                        D1O   vector of length  unique_molecules
                                      of concentration change per unit time.
				      Set by this routine.

    choice                      I0I   0 for lr_approximate_delta_concs
                                      2 for ce_approximate_delta_concs

  */
  double alt_flux_scaling;
  int success;
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
  case 42:
    /*
      Use debugging flavor with kinetic rate constants for coupledenzyme.in
    */
    success = ce_approximate_delta_concs(state, 
					 concs,
					 flux, 
					 choice);
    break;
  default :
    success = 0;
    if (lfp) {
      fprintf(lfp,"approximate_delta_concs: invalid delta_concs_choice was %d\n",choice);
      fflush(lfp);
    }
  } /* end switch (choice) */
  return(success);
}
