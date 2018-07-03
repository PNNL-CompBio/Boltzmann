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
  /*
  struct pseudoisomer_struct *pseudoisomers;
  */
  struct reaction_struct *reactions;
  struct reaction_struct *reaction;
  struct reactions_matrix_struct *rxns_matrix;
  double rxn_dg0_tf;
  double m_dg0_tf;
  int64_t *rxn_ptrs;
  int64_t *molecules_indices;
  double *coefficients;
  int64_t *matrix_text;

  double *molecule_dg0tfs;
  int64_t print_output;

  int    *dg0tfs_set;

  double coeff;

  int nrxns;
  int rxns;

  int j;
  int k;

  int use_dgzero;
  int success;

  int use_pseudoisomers;
  int computable;

  FILE* lfp;
  FILE* efp;

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
  use_dgzero               = state->use_dgzero;
  use_pseudoisomers        = state->use_pseudoisomers;
  molecule_dg0tfs          = state->molecule_dg0tfs;
  dg0tfs_set               = state->dg0tfs_set;

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
    if (use_pseudoisomers == 0) {
      /*
	User has requested to only use read in dgzero's for reaction
	delta_g0 values. If one was not specified it is an error.
	If reaction->deltag0_computed is 0, one was read in and 
	is already in reaction->deltag0;
      */
      if (reaction->deltag0_computed < 0) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"compute_reaction_dg0: You specified USE_PSEUDOISOMERS 0,"
		  " and reaction %d did not have a valid DGZERO line.\n",rxns);
	  fflush(lfp);
	}
	break;
      } 
    } else {
      /*
	If the user has specified USE_DGZERO = 1, then let a 
	read in DGZERO supercede the computed value.
      if ((use_dgzero == 1)  && (reaction->deltag0_computed == 0)) {
	do nothing here as reaction->deltag0 has been set in
	parse_reactions_file.
	}
      */
      if ((use_dgzero == 0) || (reaction->deltag0_computed < 0)) {
	/*
	  Either or use_dgzero is not set or it was set, but no
	  DGZERO line was present for this reaction, So try to
	  compute one.
	*/
	computable = 1;
	rxn_dg0_tf = (double)0.0;
	for (j=rxn_ptrs[rxns];j<rxn_ptrs[rxns+1];j++) {
	  k = molecules_indices[j];
	  coeff = coefficients[j];
	  if (dg0tfs_set[k]) {
	    m_dg0_tf = molecule_dg0tfs[k];
	    rxn_dg0_tf += (coeff * m_dg0_tf);
	  } else {
	    computable = 0;
	    break;
	  }
	} /* end for(j...) */
	if (computable) {
	  /*
	    Set the free energy units to KJ/mol.
	  */
	  reaction->unit_i = 1;
	  reaction->delta_g0 = rxn_dg0_tf;
	  reaction->deltag0_computed = 1;
	  if (print_output) {
	    if (lfp) {
	      fprintf(lfp,"Computed reaction %d delta_g0 = %le\n",rxns,rxn_dg0_tf);
	    }
	  }
	} else {
	  /*
	    This is an error condition, as use_deltag0 was 0, or
	    it was 1 but no DGZERO line existed in the reactions file for
	    this reaction.
	  */
	  success = 0;
	  if (lfp) {
	    fprintf(lfp,"compute_reaction_dg0: For reaction %d, delta_g0 was not computable.\n",rxns);
	    if (use_dgzero == 1) {
	      fprintf(lfp,"And this reaction did not have a DGZERO line in the reactions file\n");
	    }   
	  }
	  break;
	}
      } else {
	if (print_output) {
	  fprintf(lfp,"Read in reaction %d delta_g0 = %le\n",
		  rxns,reaction->delta_g0);
	}
      } /* end else if not read in */
    } /* end else use_pseudoisomers > 0 */
    reaction += 1; /* Caution address arithmetic */
  } /* end for (rxns...) */
  if (print_output) {
    if (lfp) {
      fprintf(lfp,"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
    }
  }
  return (success);
}
