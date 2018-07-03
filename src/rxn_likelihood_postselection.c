/* rxn_likelihood_postselection.c
*******************************************************************************
boltzmann

Pacific Northwest National Laboratory, Richland, WA 99352.

Copyright (c) 2010 Battelle Memorial Institute.

Publications based on work performed using the software should include 
the following citation as a reference:


Licensed under the Educational Community License, Version 2.0 (the "License"); 
you may not use this file except in compliance with the License. 
The terms and conditions of the License may be found in 
ECL-2.0_LICENSE_TERMS.TXT in the directory containing this file.
        
Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR 
CONDITIONS OF ANY KIND, either express or implied. See the License for the 
specific language governing permissions and limitations under the License.
******************************************************************************/

#include "boltzmann_structs.h"
#include "boltzmann_cvodes_headers.h"
#include "cvodes_params_struct.h"
#include "conc_to_pow.h"
#include "rxn_likelihood_postselection.h"

double rxn_likelihood_postselection(double *counts, 
		      struct state_struct *state,
		      int rxn_direction,
		      int rxn) {
  /*
    Compute the reaction quotient Q_p as the ratio of reactant counts
    to product counts with any count where the stoichiometric
    coefficient is greater than 1 is raised to that power, and multiply by
    the reaction equilibrium constant ke.
    This could also be call rxn_likelihood?

    We want to do this in a stable fashion and we want to avoid division
    by zero. So to accomplish this in a stable fashion, the expectation
    is that all product counts would have increased and thereby
    must be nonzero as they could not be negative to start with.
    The expectation is that there will only be a few (probably fewer than 4)
    reactants or products, and thus the product of their counts
    should not over, nor under flow. If this changes we might want
    to sort reactant and product concentrations and take the product of 
    successive quotients which we would expect to be well scaled.

    Called by: metropolis
    Calls      fprintf, fflush, log (intrinsic)
  */
  struct cvodes_params_struct *cvodes_params;
  struct reactions_matrix_struct *rxns_matrix;
  /*
  struct species_matrix_struct species_matrix;
  */
  double rxn_likelihood;
  double *ke;
  double *kss;
  double *rcoef;
  /*
  double *kssr;
  */
  double  count;
  double  left_counts;
  double  right_counts;
  double  eq_k;
  double  coeff;
  double  factorial;
  int64_t *rxn_ptrs;
  int64_t *molecules_indices;
  int success;
  int nrxns;

  int i;
  int j;

  int use_deq;
  int padi;

  int compute_sensitivities;
  int ode_solver_choice;

  left_counts = 1.0;
  right_counts = 1.0;

  success           = 1;
  nrxns             = (int)state->number_reactions;
  rxns_matrix       = state->reactions_matrix;
  rxn_ptrs          = rxns_matrix->rxn_ptrs;
  rcoef             = rxns_matrix->coefficients;
  molecules_indices = rxns_matrix->molecules_indices;
  ke                = state->ke;
  kss               = state->kss;
  use_deq           = state->use_deq;
  ode_solver_choice = state->ode_solver_choice;
  compute_sensitivities = state->compute_sensitivities;
  if (use_deq && compute_sensitivities && (ode_solver_choice == 1)) {
    cvodes_params = state->cvodes_params;
    ke            = cvodes_params->p;
  }
  /*
  kssr              = state->kssr;
  */
  i = rxn;
  left_counts        = 1.0;
  right_counts       = 1.0;
  eq_k = ke[rxn] * kss[rxn];
  /*
    This may change if kssr[rxn] != 1/kss[rxn]
  */
  if (rxn_direction < 0) {
    eq_k = 1.0/eq_k;
  }
  for (j=rxn_ptrs[i];j<rxn_ptrs[i+1];j++) {
    coeff = rcoef[j];
    count  = counts[molecules_indices[j]];
    if (rxn_direction < 0) {
      coeff = 0.0 -coeff;
    }
    if (coeff < 0.0) {
      coeff = 0.0 - coeff;
      factorial = -1.0;
      left_counts = left_counts * conc_to_pow(count,coeff,factorial);
      /*
      for (k=0;k<(0-coeff);k++) {
	left_counts = left_counts * (count-k);	
      } 
      */
    } else {
      /*
	NB if the coeff is == 0 then this molecule was a solvent and
	is not to be used in computing likelihoods.
      */
      if (coeff > 0.0) {
	factorial = 1.0;
	right_counts = right_counts * conc_to_pow(count,coeff,factorial);
	/*
	for (k=1;k<=coeff;k++) {
	  right_counts = right_counts * (count+k);
	} 
	*/
      }
    }
  }
  rxn_likelihood = eq_k * (left_counts/ right_counts);
  return(rxn_likelihood);
}
