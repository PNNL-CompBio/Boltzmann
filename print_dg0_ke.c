/* print_dg0_ke.c
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

#include "print_dg0_ke.h"
int print_dg0_ke(struct state_struct *state) {
  /*
    Print the reaction delta G_0's and equilibrium constants, k_e, along with
    the reaction title and stoichiometric statement.
    Called by: echo_inputs
    Calls:     fopen, fprintf, fclose (intrinsic)
  */
  struct reaction_struct *reaction;
  struct reactions_matrix_struct *rxns_matrix;
  int64_t *rxn_ptrs;
  int64_t *coefficients;
  int64_t *matrix_text;
  double *dg0s;
  double *ke;

  char *molecules_text;
  char *rxn_title_text;
  char *title;
  char *molecule;

  int success;
  int rxns;

  int nrxns;
  int np;

  int nr;
  int coeff;

  int j;
  int padi;

  FILE *dg0_ke_fp;

  success = 1;
  molecules_text   = state->molecules_text;
  rxn_title_text   = state->rxn_title_text;

  dg0_ke_fp = fopen(state->dg0ke_file,"w+");
  if (dg0_ke_fp == NULL) {
    fprintf(stderr,
	    "print_dg0_ke.c: Error could not open %s file.\n",
	    state->dg0ke_file);
    success = 0;
  }
  if (success) {
    reaction       	 = state->reactions;
    rxns_matrix    	 = state->reactions_matrix;
    rxn_ptrs       	 = rxns_matrix->rxn_ptrs;
    coefficients   	 = rxns_matrix->coefficients;
    matrix_text    	 = rxns_matrix->text;
    nrxns                = (int)state->number_reactions;
    ke                   = state->ke;
    dg0s                 = state->dg0s;

    fprintf(dg0_ke_fp,"Rxn_title\tReaction\tDelta G0\tKe\n");
    for (rxns=0;rxns < nrxns;rxns++) {
      if (reaction->title>=0) {
	title = (char *)&rxn_title_text[reaction->title];
	fprintf(dg0_ke_fp,"%s\t",title);
      }
      nr = 0;
      for (j=rxn_ptrs[rxns];j<rxn_ptrs[rxns+1];j++) {
	coeff = coefficients[j];
	if (coeff < 0) {
	  if (coeff < -1) {
	    coeff = -coeff;
	    fprintf(dg0_ke_fp,"%d ",coeff);
	  }
	  molecule = (char*)&molecules_text[matrix_text[j]];
	  fprintf(dg0_ke_fp,"%s",molecule);
	  nr += 1;
	  if (nr < reaction->num_reactants) {
	    fprintf(dg0_ke_fp," + ");
	  } else {
	    fprintf(dg0_ke_fp," => ");
	    break;
	  }
	}
      }
      np = 0;
      for (j=rxn_ptrs[rxns];j<rxn_ptrs[rxns+1];j++) {
	coeff = coefficients[j];
	if (coeff > 0) {
	  if (coeff > 1) {
	    fprintf(dg0_ke_fp,"%d ",coeff);
	  }
	  molecule = (char*)&molecules_text[matrix_text[j]];
	  fprintf(dg0_ke_fp,"%s",molecule);
	  np += 1;
	  if (np < reaction->num_products) {
	    fprintf(dg0_ke_fp," + ");
	  } else {
	    fprintf(dg0_ke_fp," ");
	    break;
	  }
	}
      }
      fprintf(dg0_ke_fp,"\t%le\t%le\n",dg0s[rxns],ke[rxns]);
      reaction += 1; /* Caution address arithmetic */
    }
    fclose(dg0_ke_fp);
  }
  return(success);
}
