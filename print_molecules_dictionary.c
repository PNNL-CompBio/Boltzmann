/* print_molecules_dictionary.c
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

#include "print_molecules_dictionary.h"
int print_molecules_dictionary(struct state_struct *state) {
  /*
    Print the  unique molecules in sorted order.
    Called by: boltzmann_init
    Calls:     fopen, fprintf, fclose (intrinsic)
  */
  struct rxn_struct *reactions;
  struct rxn_matrix_struct *rxns_matrix;
  int64_t *column_indices;
  struct istring_elem_struct *cur_molecules;
  struct istring_elem_struct *cur_cmpts;
  struct istring_elem_struct *cur_cmpt;
  char *cmpt_string;
  int nzr;
  int i;

  int oi;
  int ci;

  int success;
  int nu_molecules;

  FILE *dict_fp;
  success = 1;
  nu_molecules     = state->unique_molecules;
  cur_molecules = state->sorted_molecules;
  cur_cmpts     = state->sorted_cmpts;
  dict_fp = fopen("molecules.dict","w+");
  if (dict_fp == NULL) {
    fprintf(stderr,
	    "print_molecules_dictionary: Error could not open "
	    "molecules.dict file.\n");
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
      if (ci != -1) {
	fprintf(dict_fp,"%d %s %s\n",i,cur_molecules->string,cmpt_string);
      } else {
	fprintf(dict_fp,"%d %s\n",i,cur_molecules->string);
      }
      cur_molecules += 1; /* Caution address arithmetic. */
    }
    fclose(dict_fp);
  }
  return(success);
}
