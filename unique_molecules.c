/* unique_molecules.c
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

#include "unique_molecules.h"
int unique_molecules(struct state_struct *state) {
  /*
    Remove duplicates from the sorted_molecules list
    and set the column_indices fields in the
    reactions_matrix appropriately.
    Called by: boltzmann
    Calls:     strcmp (intrinsic)
  */
  struct rxn_struct *reactions;
  struct rxn_matrix_struct *rxns_matrix;
  int64_t *molecules_indices;
  struct istring_elem_struct *sorted_molecules;
  struct istring_elem_struct *cur_molecule;
  struct istring_elem_struct *umolecules_next;
  int nzr;
  int i;
  int success;
  int nu_molecules;
  success = 1;
  nzr            = state->number_molecules;
  sorted_molecules = state->sorted_molecules;
  rxns_matrix    = state->reactions_matrix;
  molecules_indices = rxns_matrix->molecules_indices;
  /* loop over sorted molecules. */
  nu_molecules = 0;
  molecules_indices[sorted_molecules->index] = nu_molecules;
  cur_molecule = sorted_molecules;
  sorted_molecules += 1; /* Caution address arithmetic. */
  umolecules_next  = sorted_molecules;
  for (i=1;i<nzr;i++) {
    if (strcmp(sorted_molecules->string,cur_molecule->string) != 0) {
      nu_molecules += 1;
      cur_molecule = sorted_molecules;
      umolecules_next->string = cur_molecule->string;
      umolecules_next += 1; /* Caution address arithmetic. */
    }
    molecules_indices[sorted_molecules->index] = nu_molecules;
    sorted_molecules += 1; /* Caution address arithmetic. */
  }
  state->unique_molecules = nu_molecules + 1;
  return(success);
}
