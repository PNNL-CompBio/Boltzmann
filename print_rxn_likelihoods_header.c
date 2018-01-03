/* print_rxn_likelihoods_header.c
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

#include "djb_timing_b.h"
#include "boltzmann_structs.h"

#include "print_rxn_likelihoods_header.h"
void print_rxn_likelihoods_header(struct state_struct *state) {
  /*
    Print the header lines for the reaction likelihoods output file.
    Called by: boltzmann_init .
    Calls    : fprintf,fflush.
  */
  struct rxn_struct *reactions;
  int i;
  int padi;
  if (state->rxn_lklhd_fp) {
    fprintf(state->rxn_lklhd_fp,"iter\tentropy\tdg_forward\tforward_rxn_likelihood\treverse_rxn_likelihood\n");
    fprintf(state->rxn_lklhd_fp,"iter\tentropy\tdg_forward");
    reactions                   = state->reactions;
    for (i=0;i<state->number_reactions;i++) {
      fprintf(state->rxn_lklhd_fp,"\tf_%s\tr_%s",reactions->title,reactions->title);
      reactions += 1; /* Caution address arithmetic */
    }
    fprintf(state->rxn_lklhd_fp,"\n");
    fflush(state->rxn_lklhd_fp);
  }
  return;
}
