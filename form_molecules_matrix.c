/* form_molecules_matrix.c
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

#include "form_molecules_matrix.h"
int form_molecules_matrix(struct state_struct *state,
			  struct molecules_matrix_struct *molecules_matrix,
			  int64_t *transpose_rp) {
  /*
    Transpose the reactions matrix in to the molecules matrix.
    Called by: boltzmann_init
    Calls:     
  */
  struct rxn_matrix_struct *rxn_matrix;
  /*
  struct molecules_matrix_struct *molecules_matrix;
  */
  int64_t *rxn_ptrs;
  int64_t *molecules_ptrs;
  int64_t *rcoef;
  int64_t *scoef;
  int64_t *molecules_indices;
  int64_t *rxn_indices;
  /*
  int64_t *transpose_rp;
  */
  int nu_molecules;
  int success;

  int nzr;
  int i;

  int j;
  int k;

  int m;
  int padi;
  success  = 1;
  nu_molecules = state->nunique_molecules;
  nzr        = state->number_molecules;
  rxn_matrix = state->reactions_matrix;
  /*
  molecules_matrix = state->molecules_matrix;
  */
  rxn_ptrs        = rxn_matrix->rxn_ptrs;
  molecules_indices = rxn_matrix->molecules_indices;
  rcoef           = rxn_matrix->coefficients;
  molecules_ptrs    = molecules_matrix->molecules_ptrs;
  rxn_indices     = molecules_matrix->rxn_indices;
  scoef           = molecules_matrix->coefficients;
  /*
  transpose_rp    = state->transpose_work;
  */
  for (i=0;i<nu_molecules;i++) {
    transpose_rp[i]  = 0;
  }
  for (i=0;i<nzr;i++) {
    j = molecules_indices[i];
    transpose_rp[j] += 1;
  }
  transpose_rp[nu_molecules] = 0;
  for (i=1;i< nu_molecules;i++) {
    transpose_rp[i] += transpose_rp[i-1];
  }
  molecules_ptrs[0] = 0;
  for (i=1;i<=nu_molecules;i++) {
    molecules_ptrs[i] = transpose_rp[i-1];
  }
  for (i=0;i<nu_molecules;i++) {
    transpose_rp[i] = molecules_ptrs[i];
  }
  for (i=0;i<nu_molecules;i++) {
    for (j=rxn_ptrs[i];j<rxn_ptrs[i+1];j++) {
      k = molecules_indices[j];
      m = transpose_rp[k];
      rxn_indices[m] = i;
      scoef[m]     = rcoef[j];
      transpose_rp[k] += 1;
    }
  }
  return(success);
}
