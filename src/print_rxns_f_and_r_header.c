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

#include "print_rxns_f_and_r_header.h"
void print_rxns_f_and_r_header(struct state_struct *state, char *index_str,
			       FILE *ofp) {
  /*
    Print a  header lines for reaction related 
    quantities forward and reverse.
    Called by: ode_print_lklhd_header, ode_print_kq_header, 
               ode_print_sqk_header
    Calls    : fprintf,fflush.
  */
  struct reaction_struct *reactions;
  char *rxn_title_text;
  char *title;
  int i;
  int nrxns;
  rxn_title_text = state->rxn_title_text;
  nrxns = state->number_reactions;
  if (ofp) {
    fprintf(ofp,"%s",index_str);
    reactions                   = state->reactions;
    for (i=0;i<nrxns;i++) {
      title = (char*)&rxn_title_text[reactions->title];
      fprintf(ofp,"\tf_%s\tr_%s",title,title);
      reactions += 1; /* Caution address arithmetic */
    }
    fprintf(ofp,"\n");
    fflush(ofp);
  }
  return;
}
