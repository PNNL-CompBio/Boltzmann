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

#include "boltzmann_structs.h"
#include "unique_compartments_core.h"

#include "unique_compartments.h"
int unique_compartments(struct state_struct *state) {
  /*
    Remove duplicates from the sorted_compartments list
    and set the compartment_indices fields in the
    reactions_matrix appropriately.
    Called by: rxns_init
    Calls:     unique_compartments_core
  */
  struct rxn_matrix_struct *rxns_matrix;
  int64_t *compartment_indices;
  struct compartment_struct *sorted_cmpts;
  char *compartment_text;
  int64_t sum_compartment_len;
  int64_t nunique_compartments;
  int64_t align_len;
  int64_t align_mask;

  int nzr;
  int success;

  success = 1;
  align_len     = state->align_len;
  align_mask     = state->align_mask;
  nzr            = state->number_compartments;
  sorted_cmpts   = state->sorted_cmpts;
  rxns_matrix    = state->reactions_matrix;
  compartment_text = state->compartment_text;
  compartment_indices = rxns_matrix->compartment_indices;
  
  success = unique_compartments_core(nzr,
				     sorted_cmpts,
				     compartment_text,
				     compartment_indices,
				     &nunique_compartments,
				     &sum_compartment_len,
				     align_len,
				     align_mask);
  state->nunique_compartments = nunique_compartments;
  state->sum_compartment_len  = sum_compartment_len;
  return(success);
}
