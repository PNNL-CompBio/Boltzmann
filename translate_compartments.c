/* translate_compartments.c
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

#include "translate_compartments.h"
int translate_compartments(struct state_struct *state) {
  /*
    Assign the proper compartment numbers to the 
    unsorted molecules, using the compartment_indices field
    of the rxns_matrix as set by call to unique_compartments.
    Called by: rxns_init
  */
  struct rxn_matrix_struct *rxns_matrix;
  struct molecule_struct *molecule;
  int64_t *compartment_indices;
  int nzr;
  int i;
  int success;
  int c_indx;
  success = 1;
  nzr         = state->number_molecules;
  molecule    = state->unsorted_molecules;
  rxns_matrix = state->reactions_matrix;
  compartment_indices = rxns_matrix->compartment_indices;
  /* loop over unsorted molecules. */
  for (i=0;i<nzr;i++) {
    c_indx = molecule->c_index;
    if (c_indx > 0) {
      molecule->c_index = compartment_indices[c_indx];
    } 
    molecule += 1; /* Caution Address arithmetic */
  }
  return(success);
}
