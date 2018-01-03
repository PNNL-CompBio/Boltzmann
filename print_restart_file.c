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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <float.h>
#include <signal.h>
#include <unistd.h>

#include "boltzmann_structs.h"

#include "print_restart_file.h"
int print_restart_file(struct state_struct *state) {
  /*
    Print the  restart file concentrations.
    Called by: boltzmann_run_sim
    Calls:     fopen, fprintf, fclose (intrinsic)
  */
  struct rxn_struct *reactions;
  struct rxn_matrix_struct *rxns_matrix;
  struct istring_elem_struct *cur_molecules;
  struct istring_elem_struct *cur_cmpts;
  struct istring_elem_struct *cur_cmpt;
  double *cconcs;
  int64_t *column_indices;
  char *cmpt_string;
  char vbs[8];
  char *vbsp;
  char *vbsc;
  char *vbsv;
  int nzr;
  int i;

  int oi;
  int ci;

  int success;
  int nu_molecules;

  FILE *restart_fp;
  success = 1;
  nu_molecules  = state->unique_molecules;
  cur_molecules = state->sorted_molecules;
  cur_cmpts     = state->sorted_cmpts;
  cconcs        = state->current_concentrations;
  restart_fp    = state->restart_fp;
  vbs[0]        = 'C';
  vbs[1]        = '\0';
  vbsc          = (char *)&vbs[0];
  vbsv          = (char *)&vbs[1];
  restart_fp = fopen(state->restart_file,"w+");
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
    for (i=0;i<nu_molecules;i++) {
      ci = cur_molecules->c_index;
      if (ci != oi) {
	oi = ci;
	cur_cmpt = (struct istring_elem_struct *)&(cur_cmpts[ci]);
	cmpt_string = cur_cmpt->string;
      }
      if (cur_molecules->variable) {
	vbsp = vbsv;
      } else {
	vbsp = vbsc;
      }
      if (ci != -1) {
	fprintf(restart_fp," %s:%s\t%le\t%s\n",
		cur_molecules->string,cmpt_string,cconcs[i],vbsp);
      } else {
	fprintf(restart_fp," %s\t%le\t%s\n",
		cur_molecules->string,cconcs[i],vbsp);
      }
      cur_molecules += 1; /* Caution address arithmetic. */
    }
    fclose(restart_fp);
  }
  return(success);
}
