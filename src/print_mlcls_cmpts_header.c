/* print_mlcls_cmpts_header.c
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
#include "print_mlcls_cmpts_header.h"
void print_mlcls_cmpts_header(struct state_struct *state,
			      char *index_str,
			      FILE *ofp) {
  /*
    print sorted molecule names with compartments (separated by :)
    title line.
    Called by: ode_print_grad_header,
               ode_print_concs_header,
    Calls:     fopen, fprintf, fflush
  */
  struct molecule_struct *cur_molecule;
  struct compartment_struct *cur_cmpts;
  struct compartment_struct *cur_cmpt;
  char *cmpt_string;
  char *molecule;
  char *molecules_text;
  char *compartment_text;
  char *ode_grad_file;
  int i;
  int ci;

  int nu_molecules;
  int padi;

  FILE *lfp;
  FILE *efp;
  nu_molecules     = state->nunique_molecules;
  cur_molecule     = state->sorted_molecules;
  cur_cmpts        = state->sorted_compartments;
  molecules_text   = state->molecules_text;
  compartment_text = state->compartment_text;
  ode_grad_file    = state->ode_grad_file;
  lfp              = state->lfp;

  if (ofp) {
    fprintf(ofp,"%s",index_str);
    for (i=0;i<nu_molecules;i++) {
      if ((cur_molecule->solvent == 0) || (cur_molecule->variable == 1)) {
	ci = cur_molecule->c_index;
	molecule = (char*)&molecules_text[cur_molecule->string];
	fprintf(ofp,"\t%s",molecule);
	if (ci >= 0) {
	  cur_cmpt   = (struct compartment_struct *)&(cur_cmpts[ci]);
	  cmpt_string = (char*)&compartment_text[cur_cmpt->string];
	  fprintf(ofp,":%s",cmpt_string);
	}
      }
      cur_molecule += 1; /* caution address arithmetic.*/
    }
    fprintf(ofp,"\n");
    fflush(ofp);
  }
}
