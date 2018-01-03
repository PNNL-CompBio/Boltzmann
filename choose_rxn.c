/* choose_rxn.c
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
#include "candidate_rxn.h"
#include "bndry_flux_update.h"
#include "metropolis.h"

#include "choose_rxn.h"
int choose_rxn(struct state_struct *state,
	       double *r_sum_likelihoodp) {
  /*
    Choose a reaction using Metropolis Monte Carlo Methods.

    Called by: boltzmann_run
    Calls:     candidate_rxn,
               bndry_flux_update
	       metropolis

    Arguments: 
    Name        TMF   Description
    state       G*B   Pointer to the state structure.
                      Modifies forward_rxn_likelihood and
		               reverse_rxn_likelihood fields,
			       future_counts
			       bndry_flux_counts 
			          (via the bndry_flux_update call)
		      Uses number_reactions, activities,
		           vgrng2_state

    r_sum_likelihoodp
                D*O   Address of a double which is set to the reciprocal of 
		      (1+sum of the forward and reverse reaction likelihoods),
		      the likelihood of no reaction being selected.
  */
  struct vgrng_state_struct *vgrng2_state;
  double *forward_rxn_likelihood;
  double *reverse_rxn_likelihood;
  double *future_counts;
  double scaling;

  int number_reactions;
  int number_reactions_t2;

  int accept;
  int i;

  int rxn_choice;
  int rxn_direction;

  int j;
  int use_metropolis;

  int success;
  int padi;
  
  success                = 1;
  number_reactions       = (int)state->number_reactions;
  forward_rxn_likelihood = state->forward_rxn_likelihood;
  reverse_rxn_likelihood = state->reverse_rxn_likelihood;
  future_counts          = state->future_counts;
  vgrng2_state           = state->vgrng2_state;
  use_metropolis         = (int)state->use_metropolis;
  number_reactions_t2    = number_reactions << 1;
  accept = 0;
  for (j=0;((j<number_reactions)&&(accept == 0));j++) {
    /*
      Get a trial choice. This sets the future_counts field of 
      the state structure.
    */
    rxn_choice = candidate_rxn(state,&scaling,r_sum_likelihoodp);
    if ((1.0 - *r_sum_likelihoodp) < state->epsilon) {
    	/*
    	  If no reaction other than the no-op is very likely 
    	  return a -1 for reaction choice.
    	*/
      accept = 0;
      break;
    }
    /*
      Compute the likelihood for this reaction.
      -1 for reverse, 1 for forward;
    */
    rxn_direction = 1;
    i = rxn_choice;
    if (rxn_choice >= number_reactions) {
      rxn_direction = -1;
      i = i - number_reactions;
    }
    /*
      Testing with no rejection.
    return(rxn_choice);
    */
    accept = 1;
    if (rxn_choice < number_reactions_t2) {
      /*
    	The Metropolis method used to follow here.
    	If reaction was forward or reverse 
    	( No-op is rxn_choice >= number_reactions_t2).
      */
      if (use_metropolis) {
	accept = metropolis(state,rxn_direction,i,scaling);
      } else {
	/*
	  Update the boundary fluxes.
	*/
	success = bndry_flux_update(i,rxn_direction,state);
      }
    } 
  } /* end for (j...) */
  if (accept == 0) {
    /*
      No reactions were accepted. Print error message.
    */
    fprintf(stderr,"choose_rxn: Error no likely reactions.\n");
    fflush(stderr);
    rxn_choice = -1;
  } 
  return(rxn_choice);
}

