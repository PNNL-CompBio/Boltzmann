/* ode_print_lklhd_header.c
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

#include "ode_print_lklhd_header.h"
void ode_print_lklhd_header(struct state_struct *state) {
  /*
    Print the header lines for the ode reaction likelihoods output file.
    Called by: deq_run .
    Calls    : fprintf,fflush.
  */
  struct rxn_struct *reactions;
  char *rxn_title_text;
  char *title;
  int i;
  int nrxns;
  FILE *ode_lklhd_fp;
  FILE *lfp;
  rxn_title_text = state->rxn_title_text;
  ode_lklhd_fp   = fopen(state->ode_lklhd_file,"w");
  state->ode_lklhd_fp = ode_lklhd_fp;
  nrxns = state->number_reactions;
  if (ode_lklhd_fp) {
    fprintf(ode_lklhd_fp,"Reaction Likelihoods\n");
    fprintf(ode_lklhd_fp,"time");
    reactions                   = state->reactions;
    for (i=0;i<nrxns;i++) {
      title = (char*)&rxn_title_text[reactions->title];
      fprintf(ode_lklhd_fp,"\tf_%s\tr_%s",title,title);
      reactions += 1; /* Caution address arithmetic */
    }
    fprintf(ode_lklhd_fp,"\n");
  }
  return;
}
