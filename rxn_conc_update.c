/* rxn_conc_update.c
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

#include "rxn_conc_update.h"

int rxn_conc_update(int rxn_no, int direction,
		    struct state_struct *state) {
  /*
    Compute the change in concentrations for forward reaction rxn_no.
    Called by candidate_rxn.
    Calls:
  */
  struct istring_elem_struct *sorted_molecules;
  struct istring_elem_struct *molecule;
  double *current_concs;
  double *future_concs;
  int64_t *rxn_ptrs;
  int64_t *molecules_indices;
  int64_t *coefficients;

  int     nu_molecules;
  int     i;

  int     j;
  int     k;

  int     success;
  int     padi;
			 
  struct rxn_matrix_struct *rxns_matrix;

  success           = 1;
  current_concs     = state->current_concentrations;
  future_concs      = state->future_concentrations;
  nu_molecules      = state->nunique_molecules;
  sorted_molecules  = state->sorted_molecules;
  rxns_matrix       = state->reactions_matrix;
  rxn_ptrs          = rxns_matrix->rxn_ptrs;
  molecules_indices = rxns_matrix->molecules_indices;
  coefficients      = rxns_matrix->coefficients;
  /*
    The following could be done with an memmove
  */
  if (future_concs != current_concs) {
    for (i=0;i<nu_molecules;i++) {
      future_concs[i] = current_concs[i];
    }
  }
  if (direction > 0) {
    for (j=rxn_ptrs[rxn_no];j<rxn_ptrs[rxn_no+1];j++) {
      k = molecules_indices[j];
      molecule = (struct istring_elem_struct *) &sorted_molecules[k];
      if (molecule->variable) {
        future_concs[k] += (double)coefficients[j];
      }
      if (future_concs[k] < 0.0) {
        future_concs[k] = 0.0;
      }
    } 
  } else {
    for (j=rxn_ptrs[rxn_no];j<rxn_ptrs[rxn_no+1];j++) {
      k = molecules_indices[j];
      molecule = (struct istring_elem_struct *) &sorted_molecules[k];
      if (molecule->variable) {
        future_concs[k] -= (double)coefficients[j];
      }
      if (future_concs[k] < 0.0) {
        future_concs[k] = 0.0;
      }
    }
  }
  return(success);
}
