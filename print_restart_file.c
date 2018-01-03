/* print_restart_file.c
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

#include "print_restart_file.h"
int print_restart_file(struct state_struct *state) {
  /*
    Print the  restart file concentrations.
    Called by: boltzmann_run, deq
    Calls:     fopen, fprintf, fclose (intrinsic)
  */
  struct molecule_struct *cur_molecules;
  struct compartment_struct *cur_cmpts;
  struct compartment_struct *cur_cmpt;
  double *ccounts;
  double multiplier;
  double conc;
  char   *restart_file;
  char *cmpt_string;
  char *molecule;
  char *molecules_text;
  char *compartment_text;
  char vbs[8];
  char *vbsp;
  char *vbsc;
  char *vbsv;

  int oi;
  int ci;

  int success;
  int nu_molecules;

  int i;
  int padi;

  FILE *restart_fp;
  success = 1;
  nu_molecules  = state->nunique_molecules;
  cur_molecules = state->sorted_molecules;
  cur_cmpts     = state->sorted_cmpts;
  ccounts       = state->current_counts;
  /*multiplier    = state->count_to_conc; */
  
  molecules_text = state->molecules_text;
  compartment_text = state->compartment_text;
  vbs[0]        = 'F';
  vbs[1]        = 'V';
  vbsc          = (char *)&vbs[0];
  vbsv          = (char *)&vbs[1];
  restart_file  = state->restart_file;
  restart_fp = fopen(restart_file,"w+");
  if (restart_fp == NULL) {
    fprintf(stderr,
	    "print_restart_file: Error opening restart_file %s\n",
	    state->restart_file);
    fflush(stderr);
    success = 0;
  }
  cmpt_string = NULL;
  oi          = -1;
  if (success) {
    /*
      First print a volume, and then a CONC_UNITS line.
    */
    fprintf(restart_fp,"VOLUME          %le\n",state->default_volume);
    fprintf(restart_fp,"CONC_UNITS      %le\n",state->conc_units);
    for (i=0;i<nu_molecules;i++) {
      ci = cur_molecules->c_index;
      if (ci != oi) {
	oi = ci;
	cur_cmpt   = (struct compartment_struct *)&(cur_cmpts[ci]);
	multiplier = cur_cmpt->count_to_conc;
	if (ci > 0) {
	  cmpt_string = (char*)&compartment_text[cur_cmpt->string];
	}
      }
      if (cur_molecules->variable) {
	vbsp = vbsv;
      } else {
	vbsp = vbsc;
      }
      molecule = (char*)&molecules_text[cur_molecules->string];
      conc = ccounts[i] * multiplier;
      if (ci > 0) {
	fprintf(restart_fp," %s:%s\t%le\t%c\n",
		molecule,cmpt_string,conc,vbsp[0]);
      } else {
	fprintf(restart_fp," %s\t%le\t%c\n",
		molecule,conc,vbsp[0]);
      }
      cur_molecules += 1; /* Caution address arithmetic. */
    }
    fclose(restart_fp);
  }
  return(success);
}
