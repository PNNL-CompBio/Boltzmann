/* print_free_energy_header.c
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

#include "print_free_energy_header.h"
void print_free_energy_header(struct state_struct *state) {
  /*
    Print the header lines for the free energy output file.
    Called by: boltzmann_init .
    Calls    : fprintf,fflush.
  */
  struct reaction_struct *reactions;
  char *rxn_title_text;
  char *title;
  int i;
  int free_energy_format;

  free_energy_format = (int)state->free_energy_format;
  rxn_title_text     = state->rxn_title_text;
  /*
    Print header lines for the free energy file.
  */
  if (state->free_energy_fp) {
    if (free_energy_format == 1) {
      fprintf(state->free_energy_fp,"negative_log_likelihoods\n");
    } else if (free_energy_format == 2) {
      fprintf(state->free_energy_fp,"free energy (KJ/mol)\n");
    } else if (free_energy_format == 3) {
      fprintf(state->free_energy_fp,"free energy (Kcal/mol)\n");
    }
    fprintf(state->free_energy_fp,"iter");
    reactions = state->reactions;
    for (i=0;i<(int)state->number_reactions;i++) {
      title = (char *)&rxn_title_text[reactions->title];
      fprintf(state->free_energy_fp,"\t%s",title);
      reactions += 1; /* Caution address arithmetic */
    }
    fprintf(state->free_energy_fp,"\n");
    fflush(state->free_energy_fp);
  }
  return;
}
