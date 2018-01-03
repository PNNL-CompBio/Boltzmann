/* print_likelihoods.c
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

#include "print_likelihoods.h"
void print_likelihoods(struct state_struct *state, 
		       double entropy, 
		       double dg_forward, 
		       int step) {
  /*
    print the likelihoods
    Prints out the entropy, the delta_g0, and likelihood field of the 
    state structure in a tab delimited row terminated by a newline.

    Called by boltzmann_init, boltzmann_run_sim

    Arguments:
    
    Name          TMF      Description

    state         G*I      state structure :
                           input fields are number_reactions,
					    forward_rxn_likelihood,
					    reverse_rxn_likelihood,
					    activities,
					    rxn_lklhod_fp
                           no fields of state are modified.
			   
    entropy       DSI      double precision scalar representing system
                           entropy

    dg_forward    DSI      double presicion scalar reprenting the 
                           delta_g0 for the forward reactions.

    step          ISI      step number, >= 0.
                  
    
  */
  double *forward_rxn_likelihood;
  double *reverse_rxn_likelihood;
  double *activities;

  int j;
  int number_reactions;
  FILE *rxn_lklhd_fp;

  number_reactions     	 = (int)state->number_reactions;
  rxn_lklhd_fp         	 = state->rxn_lklhd_fp;
  forward_rxn_likelihood = state->forward_rxn_likelihood;
  reverse_rxn_likelihood = state->reverse_rxn_likelihood;
  activities        	 = state->activities;

  if (rxn_lklhd_fp) {
    fprintf(rxn_lklhd_fp,"%d\t%le\t%le",step,entropy,dg_forward);
    for (j=0;j<number_reactions;j++) {
      fprintf(rxn_lklhd_fp,"\t%le",forward_rxn_likelihood[j]*activities[j]);
      fprintf(rxn_lklhd_fp,"\t%le",reverse_rxn_likelihood[j]*activities[j]);
    }
    fprintf(rxn_lklhd_fp,"\n");
  }
  return;
}

