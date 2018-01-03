/* rxn_likelihood.c
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

#include "rxn_likelihood.h"
double rxn_likelihood(double *counts, 
		      struct state_struct *state,
		      int rxn_direction,
		      int rxn) {
  /*
    Compute the reaction quotient Q_p as the ratio of reactant counts
    to product counts with any count where the stoichiometric
    coefficient is greater than 1 is raised to that power, and multiply by
    the reaction equilibrium constant ke.

    We want to do this in a stable fashion and we want to avoid division
    by zero. So to accomplish this in a stable fashion, the expectation
    is that all product counts would have increased and thereby
    must be nonzero as they could not be negative to start with.
    The expectation is that there will only be a few (probably fewer than 4)
    reactants or products, and thus the product of their counts
    should not over, nor under flow. If this changes we might want
    to sort reactant and product counts and take the product of 
    successive quotients which we would expect to be well scaled.

    Called by: rxn_likelihoods
    Calls      fprintf, fflush, log (intrinsic)

    Arguments:
     Name           TMF          Descripton  

     counts          D*I         double precision vector of length number-
                                 unique-molecules with the molecule
				 counts to be used in the reaction
				 likelihood computation.
				 
     state          G*I		 The boltzmann state structure. No fields
                                 of this structure are modified by this
				 routine. The number_reactions, and
				 reactions fields are used as inputs,
				 and within the reactions field,
				 the rxns_ptrs, coefficients, and 
				 molecules_indices fields are used.

     rxn_direction  ISI          Scalar integer. -1 for compute the likelihood
                                 of the reverse reaction, +1 for compute the
				 likelihood of the forward reaction.

     rxn            ISI          Index of reaction for which the likelihood
                                 ratio is desired.

    Return value:
     liklehood      DSO          The reaction likelihood ratio, a double 
                                 precision scalar, for the specified reaction
				 in the specified direction given the molecule
				 counts in the counts vector.

  */
  struct rxn_matrix_struct *rxns_matrix;
  struct molecule_struct  *molecules;
  struct molecule_struct  *molecule;
  
  double likelihood;
  double *ke;
  double *kss;
  /*
  double *kssr;
  */
  double  count;
  double  left_counts;
  double  right_counts;
  double  left_concs;
  double  right_concs;
  double  eq_k;
  double  *count_to_conc;
  double  volume_recip;
  int64_t m_index;
  int64_t coeff;
  int64_t *rcoef;
  int64_t *rxn_ptrs;
  int64_t *molecules_indices;
  int success;
  int nrxns;

  int j;
  int k;

  success           = 1;
  nrxns             = (int)state->number_reactions;
  rxns_matrix       = state->reactions_matrix;
  rxn_ptrs          = rxns_matrix->rxn_ptrs;
  rcoef             = rxns_matrix->coefficients;
  molecules_indices = rxns_matrix->molecules_indices;
  molecules         = state->sorted_molecules;
  count_to_conc     = state->count_to_conc;
  ke                = state->ke;
  kss               = state->kss;
  /*
  kssr              = state->kssr;
  */
  left_counts        = 1.0;
  right_counts       = 1.0;
  left_concs = 1.0;
  right_concs = 1.0;
  eq_k = ke[rxn] * kss[rxn];
  /*
    This may change if kssr[rxn] != 1/kss[rxn]
  */
  if (rxn_direction < 0) {
    eq_k = 1.0/eq_k;
  }
  /*
    Mod here we need to account for concentrations not counts, so we
    need an left sum of coefficients and a right sum of coeffiecients.
  */
  for (j=rxn_ptrs[rxn];j<rxn_ptrs[rxn+1];j++) {
    coeff = rcoef[j];
    m_index = molecules_indices[j];
    molecule = &molecules[m_index];
    volume_recip = count_to_conc[m_index];
    /*
    */
    count = counts[m_index];
    if (rxn_direction < 0) {
      coeff = -coeff;
    }
    if (coeff < 0) {
      for (k=0;k<(0-coeff);k++) {

	left_concs = left_concs * ((count-k) * volume_recip);
	left_counts = left_counts * (count-k);

      } 
    } else {
      if (coeff > 0) {
	for (k=1;k<=coeff;k++) {

	  right_concs = right_concs * (count+k) * volume_recip;
	  right_counts = right_counts * (count+k);
	} 
      }
      /*
	NB if the coeff is == 0 then this molecule was a solvent and
	is not to be used in computing likelihoods.
      */
    }
  }
  /*
    if left_concs < 0 then not enought reatants were available to run the reaction. So set left_concs to be 0.
  */
  if (left_concs < 0.0) {
    left_concs = 0.0;
  }
  if (left_counts < 0.0) {
    left_counts = 0;
  }
  /*
  likelihood = eq_k * (left_concs/right_concs);
  */
  likelihood = eq_k * (left_counts/right_counts);
  return(likelihood);
}
