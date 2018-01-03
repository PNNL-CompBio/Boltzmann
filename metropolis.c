/* metropolis.c
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

#include "rxn_likelihood_postselection.h"
#include "vgrng.h"
#include "bndry_flux_update.h"

#include "metropolis.h"
/*
  This code used to be in choose_rxn.
  needs future_counts vector, state, rxn_direction, rxn number(i),
  scaling, reverse and forward_rxn_likelihood  and activities vectors.
*/
int metropolis(struct state_struct *state,
	       int rxn_direction,
	       int rxn_number,
	       double scaling) {
  /*
    Use the metropolis method to accept or reject a reaction
    choice based on its likelihood.
    scaling is set in candidate_rxn and is 
    the global sum of the likelihoods + 1 divided by the random number range.
    
    Returns a 1 if the reaction was accepted, and the
    boundary fluxes updated, 0 otherwise.
    
    Called by: choose_rxn
    Calls:     rxn_likelihood_postselection, vgrng, bndry_flux_update
  */
  struct vgrng_state_struct *vgrng2_state;
  double *future_counts;
  double *activities;
  double *reverse_rxn_likelihood;
  double *forward_rxn_likelihood;
  double likelihood;
  double dchoice;
  int64_t choice;
  int accept;
  int success;
  future_counts          = state->future_counts;
  forward_rxn_likelihood = state->forward_rxn_likelihood;
  reverse_rxn_likelihood = state->reverse_rxn_likelihood;
  activities             = state->activities;
  vgrng2_state           = state->vgrng2_state;
  /*
    Compute the reaction likelihood for this reaction.
  */
  likelihood = rxn_likelihood_postselection(future_counts,
					    state,rxn_direction,
					    rxn_number);
    	  
  if (likelihood < 1.0) {
    /*
      An unlikely reaction but we do give a small but nonzero
      opportunity to fire, so we grab a new random number and 
      scale it in 0:vall where vall is the global sum of the 
      likelihoods + 1 and 
      we accept the reaction if this new scaled choice is 
      less than the likelihood of the reaction. 
      Should that be likelihood * activity[rxn_choice]?;
      Doug worries that this may need a separate random number
      generator to decouple the the re-accpetance from position
      in the reaction list. Scaling is determined by call
      to candidate_rxn.
    */
    choice  = vgrng(vgrng2_state);
    dchoice = ((double)choice)*scaling;

    if (dchoice < likelihood*activities[rxn_number]) {
      accept = 1;
      /*
	Update the boundary fluxes.
      */
      success = bndry_flux_update(rxn_number,rxn_direction,state);
    } else {
      /*
	Prevent this reaction from happening again.
      */
      accept = 0;
      if (rxn_direction < 0) {
	reverse_rxn_likelihood[rxn_number] = 0.0;
      } else {
	forward_rxn_likelihood[rxn_number] = 0.0;
      }
    }
  } else {
    /* 
       new likelihood was > 1 so accept. 
    */
    accept = 1;
    /*
      Update the boundary fluxes.
    */
    success = bndry_flux_update(rxn_number,rxn_direction,state);
  }
  return(accept);
}
