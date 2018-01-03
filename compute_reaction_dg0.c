/* compute_reaction_dg0.c
 ********************************************************************************
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
/*
 *  Created on: Apr 18, 2013
 *      Author:  Dennis G. Thomas
 *  Modified by Doug Baxter on Jun 1, 2013
*/
#include "boltzmann_structs.h"

#include "compute_reaction_dg0.h"

int compute_reaction_dg0(struct state_struct *state){
  /*
    Set the delta_g0 field for all of the reactions
    based on energies of formation.
    Called by: compute_standard_energies
    Calls:     fprintf, fflush
  */
  int success;
  int64_t sum;


  double rxn_dg0_tf;
  double m_dg0_tf;

  /*
  struct pseudoisomer_struct *pseudoisomers;
  */
  struct rxn_struct *reactions;
  struct rxn_struct *reaction;
  struct rxn_matrix_struct *rxns_matrix;
  int64_t *rxn_ptrs;
  int64_t *molecules_indices;
  int64_t *coefficients;
  int64_t *matrix_text;

  double *molecule_dg0tfs;
  int64_t print_output;
  int64_t mto;

  char   *molecule;
  char   *molecules_text;

  int nrxns;
  int rxns;

  int j;
  int k;

  int coeff;
  int padi;


  FILE* lfp;


  success = 1;
  
  reactions       	   = state->reactions;
  reaction                 = reactions;
  rxns_matrix    	   = state->reactions_matrix;
  print_output             = state->print_output;
  rxn_ptrs       	   = rxns_matrix->rxn_ptrs;
  molecules_indices        = rxns_matrix->molecules_indices;
  coefficients   	   = rxns_matrix->coefficients;
  matrix_text    	   = rxns_matrix->text;
  nrxns                    = (int)state->number_reactions;
  molecules_text           = state->molecules_text;
  
  molecule_dg0tfs          = state->molecule_dg0tfs;

  lfp                  = state->lfp;
  print_output         = print_output && lfp;
  if (print_output) {
    fprintf(lfp, "Output from compute_reaction_dg0.c: \n");
    fprintf(lfp,"compute_reaction_dg0: number of reactions = %i\n",nrxns);
  }

  for (rxns=0;rxns < nrxns;rxns++) {
    if (print_output) {
      fprintf(lfp,"reaction no. %i\n",rxns);
      fprintf(lfp,"Number of reactants = %i\n",reaction->num_reactants);
      fprintf(lfp,"Number of products = %i\n",reaction->num_products);
    }
    /*
      Set the free energy units to KJ/mol.
    */
    reaction->unit_i = 1;
    sum = 0;
    rxn_dg0_tf = (double)0.0;
    for (j=rxn_ptrs[rxns];j<rxn_ptrs[rxns+1];j++) {
      k = molecules_indices[j];
      coeff = coefficients[j];
      m_dg0_tf = molecule_dg0tfs[k];
      mto = matrix_text[j];
      molecule = (char*)&molecules_text[mto];
      if (m_dg0_tf != 0.0) {
	if (print_output) {
	  fprintf(lfp,"coefficient of molecule %s = %i and dg0_tf = %f kJ/mol.\n",
		  molecule,coeff,m_dg0_tf);
	}
	rxn_dg0_tf += (coeff * m_dg0_tf);
      } else {
	sum += 1;
	rxn_dg0_tf = (double)0.0;
	if (print_output) {
	  fprintf(lfp,"molecule %s reaction %d, did not have a "
		  "formation energy.\n",molecule,rxns);
	  fflush(lfp);
	}
      }
    } /* end for(j...) */
    if (sum == 0) {
      reaction->delta_g0 = rxn_dg0_tf;
      if (print_output) {
	fprintf(lfp,"Computed reaction delta_g0 = %le\n",rxn_dg0_tf);
      }
    } else {
      if (print_output) {
	fprintf(lfp,"Using input reaction delta_g0 = %le\n",
		reaction->delta_g0);
      }
    }
    reaction += 1; /* Caution address arithmetic */
  } /* end for (rxns...) */
  if (print_output) {
    fprintf(lfp,"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
  }
  return (success);
}
