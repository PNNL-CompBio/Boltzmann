/* candidate_rxn.c
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
#include "vgrng.h"
#include "binary_search_l_u_b.h"
#include "rxn_count_update.h"

#include "candidate_rxn.h"
int candidate_rxn(struct state_struct *state, double *scalingp, 
		  double *r_sum_likelihoodp) {
  /*
    Generate a candidate reaction, to be tested by
    choose_rxn, and update counts as though
    that reaction had been selected.
    Called by : choose_rxn
    Calls     : vgrng, 
                binary_search_l_u_b,
		rxn_count_update
  */
  struct vgrng_state_struct *vgrng_state;
  double *rxn_likelihood_ps;
  double *forward_rxn_likelihood;
  double *reverse_rxn_likelihood;
  double *activities;
  double r_sum_likelihood;
  double dchoice;
  double uni_multiplier;
  double vall;
  double scaling;
  int64_t choice;
  int num_rxns;
  int num_rxns_t2;

  int num_rxns_t2_p1;
  int success;

  int rxn_choice;
  int j;

  int direction;
  int padi;

  success = 1;
  rxn_likelihood_ps      = state->rxn_likelihood_ps;
  num_rxns               = (int)state->number_reactions;
  forward_rxn_likelihood = state->forward_rxn_likelihood;
  reverse_rxn_likelihood = state->reverse_rxn_likelihood;
  activities             = state->activities;
  vgrng_state            = state->vgrng_state;
  uni_multiplier         = vgrng_state->uni_multiplier;
  num_rxns_t2            = num_rxns << 1;
  num_rxns_t2_p1         = num_rxns_t2 + 1;
  /*
    Compute the partial sums of the reaction likelihoods.
  */
  rxn_likelihood_ps[0] = forward_rxn_likelihood[0]*activities[0];
  for (j=1;j<num_rxns;j++) {
    rxn_likelihood_ps[j] = rxn_likelihood_ps[j-1] + 
      (forward_rxn_likelihood[j] * activities[j]);
  }
  for(j=0;j<num_rxns;j++) {
    rxn_likelihood_ps[num_rxns+j] = rxn_likelihood_ps[num_rxns-1+j] + 
      (reverse_rxn_likelihood[j] * activities[j]);
  }
  /*
    1.0 is added to the likelihoods to account for the 
    likeilhood that the state does not change.
  */
  vall = 1.0 + rxn_likelihood_ps[num_rxns+num_rxns-1];
  if (vall > 0.0) {
    r_sum_likelihood = 1.0/vall;
  } else {
    r_sum_likelihood = 1.0;
  }
  *r_sum_likelihoodp = r_sum_likelihood;
  rxn_likelihood_ps[num_rxns_t2] = vall;
  /*
    Unimultiplier is 1.0/2^31-1
  */
  scaling = vall * uni_multiplier;
  *scalingp = scaling;
  /*
    choice is a pseudo-random integer in [0,2^31-1] (inclusive).
  */
  choice  = vgrng(vgrng_state);
  dchoice = ((double)choice)*scaling;
  /*
    Find index of smallest dg_ps entry that is >= choice.
  */
  rxn_choice = binary_search_l_u_b(rxn_likelihood_ps,dchoice,num_rxns_t2_p1);
  if (rxn_choice < num_rxns) {
    direction = 1;
    success = rxn_count_update(rxn_choice,direction,state);
  } else {
    if (rxn_choice < num_rxns_t2) {
      direction = -1;
      success = rxn_count_update(rxn_choice-num_rxns,direction,state);
    } /* else do nothing, no change */
  }
  return (rxn_choice);
}
