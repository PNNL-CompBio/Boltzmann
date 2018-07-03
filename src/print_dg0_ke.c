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
#include "recover_solvent_coefficients.h"
#include "zero_solvent_coefficients.h"
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
  struct compartment_struct *compartments;
  struct compartment_struct *compartment;
  double *coefficients;
  double *dg0s;
  double *ke;
  int64_t *rxn_ptrs;
  int64_t *matrix_text;
  int64_t *compartment_indices;

  char *molecules_text;
  char *compartment_text;
  char *rxn_title_text;
  char *title;
  char *molecule_name;
  char *compartment_name;

  double coeff;

  int success;
  int rxns;

  int nrxns;
  int np;

  int nr;
  int padi;

  int j;
  int cmpt;

  FILE *dg0_ke_fp;
  FILE *lfp;

  success = 1;
  lfp              = state->lfp;
  molecules_text   = state->molecules_text;
  rxn_title_text   = state->rxn_title_text;
  compartment_text = state->compartment_text;

  dg0_ke_fp = fopen(state->dg0ke_file,"w+");
  if (dg0_ke_fp == NULL) {
    if (lfp) {
      fprintf(lfp,
	      "print_dg0_ke.c: Error could not open %s file.\n",
	      state->dg0ke_file);
      fflush(lfp);
      success = 0;
    }
  }
  if (success) {
    /*
      Recover the solvent coefficients.
    */
    success = recover_solvent_coefficients(state);
  }
  if (success) {
    reaction       	 = state->reactions;
    rxns_matrix    	 = state->reactions_matrix;
    compartments         = state->sorted_compartments;
    rxn_ptrs       	 = rxns_matrix->rxn_ptrs;
    coefficients   	 = rxns_matrix->coefficients;
    compartment_indices  = rxns_matrix->compartment_indices;
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
	cmpt  = (int)compartment_indices[j];
	compartment = (struct compartment_struct *)&compartments[cmpt];
	compartment_name = (char*)&compartment_text[compartment->string];
	if (coeff < 0.0) {
	  if (coeff != -1.0) {
	    coeff = -coeff;
	    fprintf(dg0_ke_fp,"%le ",coeff);
	  }
	  molecule_name = (char*)&molecules_text[matrix_text[j]];
	  if (cmpt > 0) {
	    fprintf(dg0_ke_fp,"%s:%s",molecule_name,compartment_name);
	  } else {
	    fprintf(dg0_ke_fp,"%s",molecule_name);
	  }
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
	cmpt  = (int)compartment_indices[j];
	compartment = (struct compartment_struct *)&compartments[cmpt];
	compartment_name = (char*)&compartment_text[compartment->string];
	if (coeff > 0.0) {
	  if (coeff != 1.0) {
	    fprintf(dg0_ke_fp,"%le ",coeff);
	  }
	  molecule_name = (char*)&molecules_text[matrix_text[j]];
	  if (cmpt > 0) {
	    fprintf(dg0_ke_fp,"%s:%s",molecule_name,compartment_name);
	  } else {
	    fprintf(dg0_ke_fp,"%s",molecule_name);
	  }
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
    /*
      Now we need to rezero the solvent coefficents.
    */
    if (success) {
      success = zero_solvent_coefficients(state);
    }
  }
  return(success);
}
