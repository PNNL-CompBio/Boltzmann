/* rxn_likelihoods.c
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <float.h>
#include <signal.h>
#include <unistd.h>

#include "boltzmann_structs.h"

#include "rxn_likelihoods.h"
int rxn_likelihoods(double *free_energy, 
		    double *concs, 
		    double *rxn_likelihood,
		    struct state_struct *state) {
  /*
    Compute the reaction quotients Q_p as the ratio of reactant concentrations
    to product concentration with any concentration where the stoichiometric
    coefficient is greater than 1 is raised to that power, and multiply by
    the reaction equilibrium constant ke.
    This could also be call rxn_likelihood?
    For some reason in the matlab code its actually called current_free_energy,
    but that just seems wrong.

    We want to do this in a stable fashion and we want to avoid division
    by zero. So to accomplish this in a stable fashion, the expectation
    is that all product concentations would have increased and thereby
    must be nonzero as they could not be negative to start with.
    The expectation is that there will only be a few (probably fewer than 4)
    reactants or products, and thus the product of their concentrations
    should not over, nor under flow. If this changes we might want
    to sort reactant and product concentrations and take the product of 
    successive quotients which we would expect to be well scaled.

    Called by: rxn_log_likelihoods
    Calls      fprintf, fflush, log (intrinsic)
  */
  struct rxn_matrix_struct *rxns_matrix;
  /*
  struct species_matrix_struct species_matrix;
  */
  double *ke;
  double small_nonzero;
  double  conc;
  double  left_concs;
  double  right_concs;
  double  top;
  double  bot;
  double  eq_k;
  int64_t coeff;
  int64_t *rcoef;
  int64_t *rxn_ptrs;
  int64_t *molecules_indices;
  /*
  int64_t *scoef;
  int64_t *molecules_ptrs;
  int64_t *rxn_indices;
  */
  int success;
  int nrxns;

  int i;
  int j;

  int k;
  int padi;
  left_concs = 1.0;
  right_concs = 1.0;

  success           = 1;
  nrxns             = state->number_reactions;
  small_nonzero     = state->small_nonzero;
  rxns_matrix       = state->reactions_matrix;
  rxn_ptrs          = rxns_matrix->rxn_ptrs;
  rcoef             = rxns_matrix->coefficients;
  molecules_indices = rxns_matrix->molecules_indices;
  ke                = state->ke;
  left_concs        = 1.0;
  right_concs       = 1.0;
  for (i=0;i<nrxns;i++) {
    for (j=rxn_ptrs[i];j<rxn_ptrs[i+1];j++) {
      coeff = rcoef[j];
      conc  = concs[molecules_indices[j]];
      if (coeff < 0) {
	for (k=0;k<(0-coeff);k++) {
	  left_concs *= conc;	
	} 
      } else {
	for (k=0;k<coeff;k++) {
	  right_concs *= conc;	
	} 
      }
    }
    eq_k = ke[i];
    /*
    if (free_energy[i] < 0) {
      top = right_concs;
      bot = left_concs;
      if (eq_k != 0.0) {
	eq_k = 1.0/eq_k;
      } else {
	fprintf(stderr,"rxn_ratio: equilbrium constant %d was zero\n",
		i);
	fflush(stderr);
	success = 0;
	break;
      }
    } else {
      top = left_concs;
      bot = right_concs;
    }
    if (top < 1) {
      top = small_nonzero;
    }
    if (bot < 1) {
      bot = small_nonzero;
    }
    rxn_ratio[i] = eq_k * (top / bot);
    */
    if (left_concs < 1.0) {
      left_concs = 1.0;
    }
    if (right_concs < 1.0) {
      right_concs = 1.0;
    }
    rxn_likelihood[i] = eq_k * (left_concs/ right_concs);
  } /* end for(i...) */
  return(success);
}
