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

#include "boltzmann_structs.h"

#include "print_molecules_dictionary.h"
int print_molecules_dictionary(struct state_struct *state) {
  /*
    Print the  unique molecules in sorted order and the 
    header line for the concs.out file.
    Called by: boltzmann_init
    Calls:     fopen, fprintf, fclose (intrinsic)
  */
  struct rxn_struct *reactions;
  struct rxn_matrix_struct *rxns_matrix;
  int64_t *column_indices;
  struct istring_elem_struct *cur_molecules;
  struct istring_elem_struct *cur_cmpts;
  struct istring_elem_struct *cur_cmpt;
  char *compartment_text;
  char *molecules_text;
  char *cmpt_string;
  char *molecule;
  int nzr;
  int i;

  int oi;
  int ci;

  int success;
  int nu_molecules;

  FILE *dict_fp;
  FILE *conc_fp;
  success = 1;
  nu_molecules     = state->nunique_molecules;
  cur_molecules = state->sorted_molecules;
  cur_cmpts     = state->sorted_cmpts;
  conc_fp       = state->concs_out_fp;
  molecules_text = state->molecules_text;
  compartment_text = state->compartment_text;
  dict_fp = fopen("molecules.dict","w+");
  if (dict_fp == NULL) {
    fprintf(stderr,
	    "print_molecules_dictionary: Error could not open "
	    "molecules.dict file.\n");
    fflush(stderr);
    success = 0;
  }
  if (conc_fp == NULL) {
    fprintf(stderr,
	    "print_molecules_dictionary: Error null concs_out_fp\n");
    fflush(stderr);
    success = 0;
  }
  cmpt_string = NULL;
  oi          = -1;
  if (success) {
    fprintf(conc_fp,"iter");
    for (i=0;i<nu_molecules;i++) {
      ci = cur_molecules->c_index;
      molecule    = (char *)&molecules_text[cur_molecules->string];
      if (ci != oi) {
	oi = ci;
	if (ci > 0) {
	  cur_cmpt = (struct istring_elem_struct *)&(cur_cmpts[ci]);
	  cmpt_string = (char *)&compartment_text[cur_cmpt->string];
	}
      }
      if (ci > 0) {
	fprintf(dict_fp,"%d %s %s\n",i,molecule,cmpt_string);
	fprintf(conc_fp,"\t%s:%s",molecule,cmpt_string);
      } else {
	fprintf(dict_fp,"%d %s\n",i,molecule);
	fprintf(conc_fp,"\t%s",molecule);
      }
      cur_molecules += 1; /* Caution address arithmetic. */
    }
    fprintf(conc_fp,"\n");
    fclose(dict_fp);
  }
  return(success);
}
