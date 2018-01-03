/* save_likelihoods.c
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

#include "save_likelihoods.h"
void save_likelihoods(struct state_struct *state, int rxn_view_pos) {
  /*
    save the likelihoods in a transposed matrix for output to the rxns.view 
    file.

    Called by: boltzmann_run_sim
    

    Arguments:
    
    Name          TMF      Description

    state         G*B      state structure :
                           input fields are number_reactions,
			                    rxn_view_hist_lngth,
					    forward_rxn_likelihood,
					    reverse_rxn_likelihood,
					    activities,
					    
					    
                           modified fields are, rxn_view_likelihoods,
			                        rev_rxn_view_likelihoods.
			                    

			   
    rxn_view_pos  ISI      which row of the saved likelihood matrix
			   we are mpdateing.
                  
    
  */

  double *rxn_view_data;
  double *rxn_view_p;
  double *rrxn_view_data;
  double *rrxn_view_p;
  double *forward_rxn_likelihood;
  double *reverse_rxn_likelihood;
  double *activities;

  int rxn_view_hist_lngth;
  int number_reactions;
  int j;
  int padi;

  number_reactions     	 = (int)state->number_reactions;
  rxn_view_hist_lngth  	 = (int)state->rxn_view_hist_lngth;
  forward_rxn_likelihood = state->forward_rxn_likelihood;
  reverse_rxn_likelihood = state->reverse_rxn_likelihood;
  activities        	 = state->activities;

  rxn_view_data     	 = state->rxn_view_likelihoods;
  rrxn_view_data    	 = state->rev_rxn_view_likelihoods;
  rxn_view_p = (double *)&rxn_view_data[rxn_view_pos];
  rrxn_view_p = (double *)&rrxn_view_data[rxn_view_pos];
  for (j = 0; j < number_reactions;j++) {
    /*
     *rxn_view_p = forward_rxn_likelihood[j];
     *rrxn_view_p = reverse_rxn_likelihood[j];
    */
    *rxn_view_p  = forward_rxn_likelihood[j]*activities[j];
    *rrxn_view_p = reverse_rxn_likelihood[j]*activities[j];
    rxn_view_p   += rxn_view_hist_lngth; /* Caution address arithmetic here. */
    rrxn_view_p  += rxn_view_hist_lngth; /* Caution address arithmetic here. */
  }
  return;
}
