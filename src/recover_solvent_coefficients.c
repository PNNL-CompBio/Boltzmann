/* recover_solvent_coefficients.c
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

#include "recover_solvent_coefficients.h"

int recover_solvent_coefficients (struct state_struct *state) {
  /*
    Recover the coefficients in the reaction matrix corresponding to
    solvent molecules, that have been zeroed out.
    Called by print_reactions_matrix, echo_reactions_file, print_dg0_ke
    Calls:
  */
  struct reactions_matrix_struct *rxns_matrix;
  struct molecule_struct *sorted_molecules;
  struct molecule_struct *molecule;
  int64_t num_molecules;
  int64_t *rcoef;
  int64_t *scoef;
  int64_t *rxn_ptrs;
  int64_t *molecules_indices;
  int64_t nrxns;
  int64_t j;
  int64_t index;
  int64_t solvent_pos;
  int success;
  int solvent_coef_count;
  int use_bulk_water;
  int padi;
  success           = 1;
  nrxns             = state->number_reactions;
  rxns_matrix       = state->reactions_matrix;
  solvent_pos       = state->solvent_pos;
  sorted_molecules  = state->sorted_molecules;
  use_bulk_water    = state->use_bulk_water;
  rxn_ptrs          = rxns_matrix->rxn_ptrs;
  num_molecules     = rxn_ptrs[nrxns];
  rcoef             = rxns_matrix->coefficients;
  scoef             = rxns_matrix->solvent_coefficients;
  molecules_indices = rxns_matrix->molecules_indices;
  solvent_coef_count = 0;
  for (j=0;j<num_molecules;j++) {
    index = molecules_indices[j];
    if (index == solvent_pos) {
      rcoef[j] = scoef[solvent_coef_count];
      solvent_coef_count += 1;
    } else {
      molecule = (struct molecule_struct *)&sorted_molecules[index];
      if (molecule->solvent && (use_bulk_water || (molecule->variable == 0))) {
	rcoef[j] = scoef[solvent_coef_count];
	solvent_coef_count += 1;
      }
    }
  }
  return(success);
}
