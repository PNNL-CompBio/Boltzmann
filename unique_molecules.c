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

#include "boltzmann_structs.h"
#include "unique_molecules_core.h"

#include "unique_molecules.h"
int unique_molecules(struct state_struct *state) {
  /*
    Remove duplicates from the sorted_molecules list
    and set the molecules_indices fields in the
    reactions_matrix appropriately.
    Called by: rxns_init
    Calls:     strcmp (intrinsic)
  */
  struct rxn_struct *reactions;
  struct rxn_matrix_struct *rxns_matrix;
  int64_t *molecules_indices;
  int64_t *compartment_indices;
  struct molecule_struct *sorted_molecules;
  struct molecule_struct *cur_molecule;
  struct molecule_struct *umolecules_next;
  char *molecules_text;
  char *solvent_string;
  char *cstring;
  char *sstring;
  int64_t sum_molecule_len;
  int64_t nunique_molecules;
  int64_t align_len;
  int64_t align_mask;

  int success;
  int nzr;

  int solvent_pos;
  int padi;
  /*
  int i;
  int nu_molecules;

  int ni;
  int cni;
  */
  success = 1;
  nzr                 = state->number_molecules;
  sorted_molecules    = state->sorted_molecules;
  molecules_text      = state->molecules_text;
  solvent_string      = state->solvent_string;
  rxns_matrix         = state->reactions_matrix;
  molecules_indices   = rxns_matrix->molecules_indices;
  align_len           = state->align_len;
  align_mask          = state->align_mask;
  success = unique_molecules_core(nzr,
				  sorted_molecules,
				  molecules_text,
				  solvent_string,
				  molecules_indices,
				  &nunique_molecules,
				  &sum_molecule_len,
				  &solvent_pos,
				  align_len,
				  align_mask);
  state->nunique_molecules = nunique_molecules;
  state->sum_molecule_len  = sum_molecule_len;
  state->solvent_pos       = (int64_t)solvent_pos;
  return(success);
}
