/* unique_compartments.c
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

#include "unique_compartments.h"
int unique_compartments(struct state_struct *state) {
  /*
    Remove duplicates from the sorted_compartments list
    and set the column_indices fields in the
    reactions_matrix appropriately.
    Called by: boltzmann_init
    Calls:     strcmp (intrinsic)
  */
  struct rxn_struct *reactions;
  struct rxn_matrix_struct *rxns_matrix;
  int64_t *compartment_indices;
  struct istring_elem_struct *sorted_cmpts;
  struct istring_elem_struct *cur_cmpt;
  struct istring_elem_struct *ucmpts_next;
  int nzr;
  int i;
  int success;
  int nu_cmpts;
  success = 1;
  nzr            = state->number_compartments;
  sorted_cmpts   = state->sorted_cmpts;
  rxns_matrix    = state->reactions_matrix;
  compartment_indices = rxns_matrix->compartment_indices;
  /* loop over sorted compartments. */
  nu_cmpts = 0;
  if (nzr > 1) {
    if (sorted_cmpts->string) {
      compartment_indices[sorted_cmpts->c_index] = nu_cmpts;
    } else {
      compartment_indices[sorted_cmpts->c_index] = -1;
    }
    cur_cmpt = sorted_cmpts;
    sorted_cmpts += 1; /* Caution address arithmetic. */
    ucmpts_next  = sorted_cmpts;
  }
  for (i=1;i<nzr;i++) {
    if (sorted_cmpts->string) {
      if ((cur_cmpt->string == NULL) ||
          (strcmp(sorted_cmpts->string,cur_cmpt->string) != 0)) {
	nu_cmpts += 1;
	cur_cmpt = sorted_cmpts;
	ucmpts_next->string = cur_cmpt->string;
	ucmpts_next += 1; /* Caution address arithmetic. */
      }
      compartment_indices[sorted_cmpts->c_index] = nu_cmpts;
    } else {
      compartment_indices[sorted_cmpts->c_index] = -1;
    }

    sorted_cmpts += 1; /* Caution address arithmetic. */
  }
  state->unique_compartments = nu_cmpts + 1;
  return(success);
}
