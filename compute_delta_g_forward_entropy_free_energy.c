/* compute_delta_g_forward_entropy_free_energy.c
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

#include "update_regulations.h"

#include "compute_delta_g_forward_entropy_free_energy.h"
int compute_delta_g_forward_entropy_free_energy(struct state_struct *state,
						double *dg_forward_p,
						double *entropy_p,
						int64_t step) {
  /*
    Compute the delta G for the forward reactions, dg_forward,
    the system entropy, and the free_energy field of state 
    (1 value per reaction).

    Called by: boltzmann_run

    Arguments         TMF          Description
    state             G*B          Pointer to the global state structure.
                                   Fields used as input only are m_rt = -RT,
				   activities,
				   forward_rxn_likelihood,
				   forward_rxn_log_likelihood
				   number_reactions.

				   Field modified is free_energy vector.

    dg_forward_p      D*O          The contents of this address are set
                                   to the total delta g of the 
				   forward reactions (sum of the
				   log likelihood ratios)

    entropy           D*O          The contents of this address are set			                           to the total entropy = negative sum of the 
                                   products of the normalized likelihoods and
				   their logs.
                                   
    step              JSI          eight byte integer recording step for 
                                   error reporting.
                                
  */
  double *forward_rxn_likelihood;
  double *forward_rxn_log_likelihood;
  /*
  double *current_counts;
  */
  double *free_energy;
  double *activities;
  double dg_forward;
  double entropy;
  double m_rt;
  double sum_likelihood;
  double r_sum_likelihood;
  double scaled_likelihood;

  int    j;
  int    number_reactions;
  int    success;
  int    use_regulation;
  /*
    Input fields.
  */
  success                    = 1;
  m_rt                       = state->m_rt;
  number_reactions           = (int)state->number_reactions;
  /*
  current_counts             = state->current_counts;
  */
  activities        	     = state->activities;
  forward_rxn_log_likelihood = state->forward_rxn_log_likelihood_ratio;
  forward_rxn_likelihood     = state->forward_rxn_likelihood;
  /*
    Output fields.
  */
  free_energy                = state->free_energy;
  use_regulation             = state->use_regulation;
  /*
    Recompute the forward reaction log likelihoods and 
    store the likelihoods in forward_rxn_likelihod 
    and their logs in forward_rxn_log_likelihood field
    of state. Now done before this call in boltzmann_run.
    success = rxn_log_likelihoods(current_counts,
                                  forward_rxn_likelihood,
				  forward_rxn_log_likelihood,
				  state,
				  forward);
  */
  if (use_regulation) {
    update_regulations(state);
  }
  dg_forward = 0.0;
  entropy = 0.0;
  if (success) {
    sum_likelihood = 0.0;
    for (j=0;j<number_reactions;j++) {
      free_energy[j] =  m_rt * forward_rxn_log_likelihood[j];
      dg_forward     += forward_rxn_log_likelihood[j];
      sum_likelihood += activities[j]*forward_rxn_likelihood[j];
    }
    dg_forward *= m_rt;
    if (sum_likelihood <= 0.0) {
      fprintf(stderr,"boltzmann_run: Error, nonpositivity sum_likelihood = %le in recording loop iteration %lld\n",sum_likelihood,step);
      success = 0;
    }
    r_sum_likelihood = 0.0;
    if (success) {
      r_sum_likelihood = 1.0/sum_likelihood;
      for (j=0;j<number_reactions;j++) {
	scaled_likelihood = activities[j]*forward_rxn_likelihood[j] * r_sum_likelihood;
	if (scaled_likelihood > 0) {
	  entropy -= scaled_likelihood * log(scaled_likelihood);
	}
      }
    }
  }
  *entropy_p    = entropy;
  *dg_forward_p = dg_forward;
  return(success);
}
