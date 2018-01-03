/* print_reactions_matrix.c
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

#include "print_reactions_matrix.h"
int print_reactions_matrix(struct state_struct *state) {
  /*
    Print the reaction likelihoods per reaction to a file along with
    the reaction title and stoichiometric statement.
    Called by: boltzmann_init
    Calls:     fopen, fprintf, fclose (intrinsic)
  */
  struct rxn_struct *reaction;
  struct rxn_matrix_struct *rxns_matrix;
  struct molecule_struct *molecules;
  struct compartment_struct *compartments;
  struct compartment_struct *cur_cmpt;
  int64_t *rxn_ptrs;
  int64_t *molecules_indices;
  int64_t *coefficients;
  int64_t *matrix_text;
  double *rxn_view_data;
  double *rev_rxn_view_data;
  double *activities;
  double *no_op_likelihood;
  int    *rxn_fire;

  char *molecules_text;
  char *compartment_text;
  char *rxn_title_text;
  char *title;
  char *molecule;
  char *cmpt_string;

  int  *mat_row;

  int success;
  int rxns;

  int np;
  int nr;
  
  int coeff;
  int j;

  int nrxns;
  int nmols;

  int ci;
  int padi;

  FILE *rxn_mat_fp;
  FILE *lfp;

  success = 1;
  molecules_text   = state->molecules_text;
  compartment_text = state->compartment_text;
  rxn_title_text   = state->rxn_title_text;

  /*
    we don'y have a rxn_mat_file filed in state yet but if we did 
    it would look like this:

  rxn_mat_fp = fopen(state->rxn_mat_file,"w+");
  if (rxn_mat_fp == NULL) {
    fprintf(stderr,
           "print_reactions_matrix: Error could not open %s file.\n",
  	   state->rxn_matrix_file);
    success = 0;
  }
  */
  rxn_mat_fp = fopen("rxns.mat","w+");
  if (rxn_mat_fp == NULL) {
    fprintf(stderr,
	    "print_reactions_matrix: Error could not open rxns.mat file.\n");
    success = 0;
  }
  if (success) {
    reaction       	 = state->reactions;
    rxns_matrix    	 = state->reactions_matrix;
    activities           = state->activities;
    rxn_ptrs       	 = rxns_matrix->rxn_ptrs;
    molecules_indices    = rxns_matrix->molecules_indices;
    coefficients   	 = rxns_matrix->coefficients;
    matrix_text    	 = rxns_matrix->text;
    nrxns                = (int)state->number_reactions;
    nmols                = (int)state->nunique_molecules;
    molecules            = state->sorted_molecules;
    compartments         = state->sorted_cmpts;
    mat_row              = state->rxn_mat_row;

    /*
      Print header line.
    */

    fprintf(rxn_mat_fp,"reaction title\tforward reaction");
    for(j=0;j<nmols;j++) {
      ci = molecules->c_index;
      molecule    = (char *)&molecules_text[molecules->string];
      if (ci > 0) {
	cur_cmpt = (struct compartment_struct *)&(compartments[ci]);
	cmpt_string = (char *)&compartment_text[cur_cmpt->string];
	fprintf(rxn_mat_fp,"\t%s:%s",molecule,cmpt_string);
      } else {
	fprintf(rxn_mat_fp,"\t%s",molecule);
      }
      molecules += 1; /* Caution address arithmetic. */
    }
    fprintf(rxn_mat_fp,"\n");

    for (rxns=0;rxns < nrxns;rxns++) {
      /*
	Initialize the matrix_row to be all zeros;
      */
      for (j=0;j<nmols;j++) {
	mat_row[j] = 0;
      }
      if (reaction->title>=0) {
	title = (char *)&rxn_title_text[reaction->title];
	fprintf(rxn_mat_fp,"%s\t",title);
      } else {
	fprintf(rxn_mat_fp," \t");
      }
      nr = 0;
      for (j=rxn_ptrs[rxns];j<rxn_ptrs[rxns+1];j++) {
	coeff = coefficients[j];
	if (coeff < 0) {
	  mat_row[molecules_indices[j]] = coeff;
	  if (coeff < -1) {
	    coeff = -coeff;
	    fprintf(rxn_mat_fp,"%d ",coeff);
	  }
	  molecule = (char*)&molecules_text[matrix_text[j]];
	  fprintf(rxn_mat_fp,"%s",molecule);
	  nr += 1;
	  if (nr < reaction->num_reactants) {
	    fprintf(rxn_mat_fp," + ");
	  } else {
	    fprintf(rxn_mat_fp," => ");
	    break;
	  }
	}
      }
      np = 0;
      for (j=rxn_ptrs[rxns];j<rxn_ptrs[rxns+1];j++) {
	coeff = coefficients[j];
	if (coeff > 0) {
	  mat_row[molecules_indices[j]] = coeff;
	  if (coeff > 1) {
	    fprintf(rxn_mat_fp,"%d ",coeff);
	  }
	  molecule = (char*)&molecules_text[matrix_text[j]];
	  fprintf(rxn_mat_fp,"%s",molecule);
	  np += 1;
	  if (np < reaction->num_products) {
	    fprintf(rxn_mat_fp," + ");
	  } else {
	    break;
	  }
	}
      }
      /*
	Now print out the coefficients (including 0's for the matrix row)
      */
      for (j=0;j<nmols;j++) {
	fprintf(rxn_mat_fp,"\t%d",mat_row[j]);
      }
      fprintf(rxn_mat_fp,"\n");
      reaction += 1; /* Caution address arithmetic here.*/
    }
    fclose(rxn_mat_fp);
  }
  return(success);
}
