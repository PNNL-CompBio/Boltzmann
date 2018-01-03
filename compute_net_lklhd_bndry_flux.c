/* compute_net_lklhd_bndry_flux.c
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

#include "compute_net_lklhd_bndry_flux.h"
void compute_net_lklhd_bndry_flux(struct state_struct *state,
				  double *net_likelihood, 
				  double *net_lklhd_bndry_flux) {
  /*
    Compute the net likehood boundary fluxes for fixed concentration
    species.
    
    Called by: ode23tb
    Calls:     

    Arguments:
    Name           TMF       Descripton
    state          G*B       state structure. 
                             Sets the net_lklhd_bndry_flux field
			     from the molecules matrix and the
			     net_likelilhood vector.
  */
  struct  molecules_matrix_struct *molecules_matrix;
  struct  molecule_struct *sorted_molecules;
  struct  molecule_struct *molecule;
  double  bndry_flux;
  int64_t *molecules_ptrs;
  int64_t *rxn_indices;
  int64_t *coefficients;
  int64_t coef;
  int64_t rxn;
  int j;
  int i;
  int unique_molecules;
  int success;

  success           = 1;
  sorted_molecules  = state->sorted_molecules;
  unique_molecules  = state->nunique_molecules;
  molecules_matrix  = state->molecules_matrix;
  molecules_ptrs    = molecules_matrix->molecules_ptrs;
  rxn_indices       = molecules_matrix->rxn_indices;
  coefficients      = molecules_matrix->coefficients;

  molecule    = sorted_molecules;
  for (i=0;i<unique_molecules;i++) {
    if (molecule->variable == 0) {
      /*
	Loop over reactions in which the molecule is involved
	summing the net_likelihoods of those where it is produced.
	(also implies subtracting the net likelihoods where it is consumed).
      */
      bndry_flux = 0.0;
      for (j=molecules_ptrs[i];j<molecules_ptrs[i+1];j++) {
	coef = coefficients[j];
	rxn  = rxn_indices[j];
	if (coef > 0) {
	  bndry_flux += net_likelihood[rxn];
	} 
	else {
	  if (coef < 0) {
	    bndry_flux -= net_likelihood[rxn];
	  }
        }
      }
      net_lklhd_bndry_flux[i] = bndry_flux;
    }
    molecule += 1; /* Caution address arithmetic.*/
  }
}
