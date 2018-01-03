/* rxn_log_likelihoods.c
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
#include "rxn_likelihoods.h"

#include "rxn_log_likelihoods.h"
int rxn_log_likelihoods(double *concs, 
			double *rxn_likelihood_values,
			double *log_rxn_likelihood,
			struct state_struct *state,
			int rxn_direction) {
  /*
    Compute the logs of the equilibrium constants and reaction quotients Q_p 
    ( = the ratio of reactant concentrations to product concentration with 
        any concentration where the stoichiometric coefficient is greater 
	than 1 is raised to that power).
    In the process we may want the rxn_likelihood as well though frequently as
    when called from update_rxn_log_likelihoods we only want the logs, so 
    rxn_likelihood_values and log_rxn_likelihood are allowed to be the same 
    vector in that case.

    We want to do this in a stable fashion and we want to avoid division
    by zero. So to accomplish this in a stable fashion, the expectation
    is that all product concentations would have increased and thereby
    must be nonzero as they could not be negative to start with.
    The expectation is that there will only be a few (probably fewer than 4)
    reactants or products, and thus the product of their concentrations
    should not over, nor under flow. If this changes we might want
    to sort reactant and product concentrations and take the product of 
    successive quotients which we would expect to be well scaled.

    Called by: update_rxn_log_likelihoods
    Calls:     rxn_likelihoods
               fprintf, fflush, log (intrinsic)
  */
  int success;
  int nrxns;
  int i;
  int padi;

  nrxns         = (int)state->number_reactions;
  success       = rxn_likelihoods(concs,rxn_likelihood_values,state,
				  rxn_direction);
  if (success) {
    for (i=0;i<nrxns;i++) {
      log_rxn_likelihood[i] = log(rxn_likelihood_values[i]);
    }
  }
  return(success);
}
