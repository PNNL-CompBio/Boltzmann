/* ode_print_concs_header.c
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
#include "ode_print_concs_header.h"
void ode_print_concs_header(struct state_struct *state) {
  /*
    Open ode_concs_file and
    Print concentration header ("Time" followed by sorted molecule names and
    compartments)
    Called by: deq_run;
    Calls:     fopen, fprintf, fflush
  */
  struct molecule_struct *cur_molecule;
  struct compartment_struct *cur_cmpts;
  struct compartment_struct *cur_cmpt;
  char *cmpt_string;
  char *molecule;
  char *molecules_text;
  char *compartment_text;
  char *ode_concs_file;
  int i;
  int ci;
  int nu_molecules;
  int padi;
  FILE *ode_concs_fp;
  FILE *lfp;
  nu_molecules     = state->nunique_molecules;
  cur_molecule     = state->sorted_molecules;
  cur_cmpts        = state->sorted_compartments;
  molecules_text   = state->molecules_text;
  compartment_text = state->compartment_text;
  ode_concs_file   = state->ode_concs_file;
  lfp              = state->lfp;

  ode_concs_fp = fopen(ode_concs_file,"w+");
  state->ode_concs_fp  = ode_concs_fp;
  if (ode_concs_fp == NULL) {
    if (lfp) {
      fprintf(lfp,"ode_print_concs_header: Error could not open %s\n",
	      ode_concs_file);
      fflush(lfp);
    }
  } else {
    fprintf(ode_concs_fp,"Time");
    for (i=0;i<nu_molecules;i++) {
      if ((cur_molecule->solvent == 0) || (cur_molecule->variable == 1)) {
      	ci = cur_molecule->c_index;
      	molecule = (char*)&molecules_text[cur_molecule->string];
      	fprintf(ode_concs_fp,"\t%s",molecule);
      	if (ci > 0) {
      		cur_cmpt   = (struct compartment_struct *)&(cur_cmpts[ci]);
      		cmpt_string = (char*)&compartment_text[cur_cmpt->string];
      		fprintf(ode_concs_fp,":%s",cmpt_string);
      	} 
      }
      cur_molecule += 1; /* caution address arithmetic.*/
    }
    fprintf(ode_concs_fp,"\n");
    fflush(ode_concs_fp);
  }
}
