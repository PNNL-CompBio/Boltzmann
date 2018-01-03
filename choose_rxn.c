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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <float.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>

#include "djb_timing_b.h"
#include "boltzmann_structs.h"
#include "vgrng.h"
#include "candidate_rxn.h"
#include "rxn_likelihood.h"

#include "choose_rxn.h"
int choose_rxn(struct state_struct *state) {
  /*
    Choose a reaction using Metropolis Monte Carlo Methods.
    Called by: boltzmann_run_sim
    Calls:     vgrng,
               candidate_rxn,
	       rxn_likelihood
    
  */
  struct vgrng_state_struct *vgrng2_state;
  double *rxn_likelihood_ps;
  double *forward_rxn_likelihood;
  double *reverse_rxn_likelihood;
  double *activities;
  double *activities_save;
  double *concs;
  double *ke;
  double dchoice;
  double scaling;
  double likelihood;
  int64_t choice;

  int num_rxns;
  int num_rxns_t2;

  int num_rxns_p1;
  int success;

  int accept;
  int i;

  int rxn_choice;
  int rxn_direction;

  int j;
  int k;
  
  int not_saved;
  int padi;

  success = 1;
  rxn_likelihood_ps      = state->rxn_likelihood_ps;
  num_rxns               = state->number_reactions;
  forward_rxn_likelihood = state->forward_rxn_likelihood;
  reverse_rxn_likelihood = state->reverse_rxn_likelihood;
  activities             = state->activities;
  activities_save        = state->activities_save;
  concs                  = state->future_concentrations;
  vgrng2_state            = state->vgrng2_state;
  num_rxns_t2            = num_rxns << 1;
  /*
    Save a copy of the activities to reset to.
  */
  accept = 0;
  not_saved  = 1;
  for (j=0;((j<num_rxns)&&(accept == 0));j++) {
    /*
      Get a trial choice. This sets the future concentrations.
    */
    rxn_choice             = candidate_rxn(state,&scaling);
    /*
      Compute the likelihood for this reaction.
      -1 for reverse, 1 for forward;
    */
    rxn_direction = 1;
    i = rxn_choice;
    if (rxn_choice >= num_rxns) {
      rxn_direction = -1;
      i = i - num_rxns;
    }
    /*
      Testing with no rejection.
    return(rxn_choice);
    */
    /*
      If reaction was forward or reverse ( No-op is rxn_choice >= num_rxns_t2).
    */
    if (rxn_choice < num_rxns_t2) {
      /*
	Compute the reaction likelihood for this reaction.
      */
      likelihood = rxn_likelihood(concs,state,rxn_direction,i);
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
	  in the reaction list.
	*/
	choice  = vgrng(vgrng2_state);
	dchoice = ((double)choice)*scaling;
	if (dchoice < likelihood*activities[i]) {
	  accept = 1;
	} else {
	  /*
	    Record existing activity levels for restoration later.
	  */
	  if (not_saved) {
	    for (k=0;k<num_rxns;k++) {
	      activities_save[k] = activities[k];
	    }
	    not_saved = 0;
	  }
	  /*
	    Prevent this reaction from happening again.
	  */
	  if (rxn_direction < 0) {
	    reverse_rxn_likelihood[i] = 0.0;
	  } else {
	    forward_rxn_likelihood[i] = 0.0;
	  }
	}
      } else {
	/* 
	  new likelihood was > 1 so accept. 
	*/
	accept = 1;
      }
    } else {
      /*
	Reaction selected was a No-op
      */
      accept = 1;
    }
  } /* end for (j...) */
  /*
    Recover the original activities if necessary.
  */
  if (not_saved == 0) {
    for (k=0;k<num_rxns;k++) {
      activities[k] = activities_save[k];
    }
  }
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

