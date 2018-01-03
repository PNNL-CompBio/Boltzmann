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
#include "rxn_likelihood.h"

#include "rxn_likelihoods.h"
int rxn_likelihoods(double *concs, 
		    double *rxn_likelihood_values,
		    struct state_struct *state,
		    int rxn_direction) {
  /*
    Compute the reaction quotients Q_p as the ratio of reactant concentrations
    to product concentration with any concentration where the stoichiometric
    coefficient is greater than 1 is raised to that power, and multiply by
    the reaction equilibrium constant ke.
    This could also call rxn_likelihood? YES.

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
    Calls      rxn_likelihood
    Arguments:
     Name           TMF          Descripton  

     concs          D*I          double precision vector of length number-
                                 unique-molecules with the molecule
				 concentrations to be use in the reaction
				 likelihood computation.
				 
     rxn_likelihood_values       double precision vector of length 
                    D*O          number_reactions with the reaction likelihood 
		                 ratios for all of the reactions in the
				 direction specified by rxn_direction.
		                 
 				 

     state          G*I		 The boltzmann state structure. No fields
                                 of this structure are modified by this
				 routine. The number_reactions, and
				 reactions fields are used as inputs,
				 and within the reactions field,
				 the rxns_ptrs, coefficients, and 
				 molecules_indices fields are used.

     rxn_direction  ISI          Scalar integer. -1 for compute the likelihood
                                 of the reverse reaction, +1 for compute the
				 likelihood of the forward reaction.

    Return value:
     success        ISO          Always 1 for this function.
  */
  int success;
  int nrxns;

  int i;
  int padi;

  success           = 1;
  nrxns             = (int)state->number_reactions;
  for (i=0;i<nrxns;i++) {
    rxn_likelihood_values[i] = rxn_likelihood(concs,state,rxn_direction,i);
  } /* end for(i...) */
  return(success);
}
