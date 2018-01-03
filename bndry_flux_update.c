/* bndry_flux_update.c
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

#include "bndry_flux_update.h"

int bndry_flux_update(int rxn_no, int direction,
		      struct state_struct *state) {
  /*
    Compute the change in boundary fluxes for reaction rxn_no,
    in the forward direction if direction > 0, otherwise in
    the reverse direction.
    Called by choose_rxn.
    Calls:
  */
  struct molecule_struct *sorted_molecules;
  struct molecule_struct *molecule;
  struct reactions_matrix_struct *rxns_matrix;
  double *bndry_flux_counts;
  int64_t *rxn_ptrs;
  int64_t *molecules_indices;
  int64_t *coefficients;

  int     nu_molecules;
  int     success;

  int     j;
  int     k;
			 
  success           = 1;
  bndry_flux_counts = state->bndry_flux_counts;
  nu_molecules      = state->nunique_molecules;
  sorted_molecules  = state->sorted_molecules;
  rxns_matrix       = state->reactions_matrix;
  rxn_ptrs          = rxns_matrix->rxn_ptrs;
  molecules_indices = rxns_matrix->molecules_indices;
  coefficients      = rxns_matrix->coefficients;
  /*
    The following could be done with an memmove
  */
  if (direction > 0) {
    for (j=rxn_ptrs[rxn_no];j<rxn_ptrs[rxn_no+1];j++) {
      k = molecules_indices[j];
      molecule = (struct molecule_struct *) &sorted_molecules[k];
      if (molecule->variable == 0) {
        bndry_flux_counts[k] += (double)coefficients[j];
      }
    } 
  } else {
    for (j=rxn_ptrs[rxn_no];j<rxn_ptrs[rxn_no+1];j++) {
      k = molecules_indices[j];
      molecule = (struct molecule_struct *) &sorted_molecules[k];
      if (molecule->variable == 0) {
        bndry_flux_counts[k] -= (double)coefficients[j];
      }
    }
  }
  return(success);
}
