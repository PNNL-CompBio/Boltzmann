/* update_rxn_log_likelihoods.c
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
#include "rxn_log_likelihoods.h"

#include "update_rxn_log_likelihoods.h"
int update_rxn_log_likelihoods(struct state_struct *state) {
  /*
    Update the forward_rxn_likelihood, forward_rxn_log_likelihood_ratio,
               reverse_rxn_likelihood, and reverse_rxn_log_likelihood_ratio
	       fields of the state vector, based on the current_counts
	       field.

    Called by: boltzmann_run
    Calls      rxn_log_likelihoods

    Arguments:
    Name        TMF       Description
    state       G*B       pointer to the state structure.
                          Input field is current_counts.
			  Modified fields are:
			  forward_rxn_likelihood
			  forward_rxn_log_likelihood_ratio
			  reverse_rxn_likelihood
			  reverse_rxn_log_likelihood_ratio
  */
  double *current_counts;
  double *forward_rxn_log_likelihood_ratio;
  double *reverse_rxn_log_likelihood_ratio;
  double *forward_rxn_likelihood;
  double *reverse_rxn_likelihood;
  int success;
  int forward;

  int reverse;
  int padi;
  
  success       		   = 1;
  forward       		   = 1;
  reverse       		   = -1;
  current_counts                   = state->current_counts;
  forward_rxn_log_likelihood_ratio = state->forward_rxn_log_likelihood_ratio;
  reverse_rxn_log_likelihood_ratio = state->reverse_rxn_log_likelihood_ratio;
  forward_rxn_likelihood           = state->forward_rxn_likelihood;
  reverse_rxn_likelihood           = state->reverse_rxn_likelihood;
  /*
  */
  rxn_log_likelihoods(current_counts,forward_rxn_likelihood,
		      forward_rxn_log_likelihood_ratio,state,forward);
  rxn_log_likelihoods(current_counts,reverse_rxn_likelihood,
		      reverse_rxn_log_likelihood_ratio,state,reverse);
  return(success);
}
