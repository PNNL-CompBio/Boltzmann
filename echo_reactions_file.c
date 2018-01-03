/* echo_reactions_file.c
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <float.h>
#include <signal.h>
#include <unistd.h>

#include "boltzmann_structs.h"

#include "echo_reactions_file.h"
int echo_reactions_file(struct state_struct *state) {
  /*
    Echo the reactions file to a rxns.echo file.
    Called by: boltzmann
  */
  struct rxn_struct *reactions;
  struct rxn_struct *reaction;
  struct rxn_matrix_struct *rxns_matrix;
  int64_t *rxn_ptrs;
  int64_t *molecules_indices;
  int64_t *coefficients;
  char **matrix_text;

  int success;
  int rxns;

  int ns;
  int padding;

  int np;
  int nr;
  
  int coeff;
  int j;

  FILE *rxn_echo_fp;
  FILE *lfp;

  char tab;
  char padc[7];

  tab  = (char)9;
  success = 1;

  rxn_echo_fp = fopen("rxns.echo","w+");
  if (rxn_echo_fp == NULL) {
    fprintf(stderr,
	    "echo_reactions_file: Error could not open rxns.echo file.\n");
    success = 0;
  }
  if (success) {
    reaction       = state->reactions;
    rxns_matrix    = state->reactions_matrix;
    rxn_ptrs       = rxns_matrix->rxn_ptrs;
    molecules_indices = rxns_matrix->molecules_indices;
    coefficients   = rxns_matrix->coefficients;
    matrix_text    = rxns_matrix->text;
    for (rxns=0;rxns < state->number_reactions;rxns++) {
      if (reaction->title) {
	fprintf(rxn_echo_fp,"REACTION\t%s\n",reaction->title);
      }
      if (reaction->pathway) {
	fprintf(rxn_echo_fp,"PATHWAY\t%s\n",reaction->pathway);
      }
      if (reaction->lcompartment) {
	if (reaction->rcompartment) {
	  fprintf(rxn_echo_fp,"LEFT_COMPARTMENT\t%s\n",reaction->lcompartment);
	  fprintf(rxn_echo_fp,
		  "RIGHT_COMPARTMENT\t%s\n",reaction->rcompartment);
	} else {
	  fprintf(rxn_echo_fp,"COMPARTMENT\t%s\n",reaction->lcompartment);
	}
      }
      fprintf(rxn_echo_fp,"LEFT\t");
      nr = 0;
      for (j=rxn_ptrs[rxns];j<rxn_ptrs[rxns+1];j++) {
	coeff = coefficients[j];
	if (coeff < 0) {
	  if (coeff < -1) {
	    coeff = -coeff;
	    fprintf(rxn_echo_fp,"%d ",coeff);
	  }
	  fprintf(rxn_echo_fp,"%s",matrix_text[j]);
	  nr += 1;
	  if (nr < reaction->num_reactants) {
	    fprintf(rxn_echo_fp," + ");
	  } else {
	    fprintf(rxn_echo_fp,"\n");
	    break;
	  }
	}
      }
      fprintf(rxn_echo_fp,"RIGHT\t");
      np = 0;
      for (j=rxn_ptrs[rxns];j<rxn_ptrs[rxns+1];j++) {
	coeff = coefficients[j];
	if (coeff > 0) {
	  if (coeff > 1) {
	    fprintf(rxn_echo_fp,"%d ",coeff);
	  }
	  fprintf(rxn_echo_fp,"%s",matrix_text[j]);
	  np += 1;
	  if (np < reaction->num_products) {
	    fprintf(rxn_echo_fp," + ");
	  } else {
	    fprintf(rxn_echo_fp,"\n");
	    break;
	  }
	}
      }
      fprintf(rxn_echo_fp,"DGZERO\t%le\n",reaction->delta_g0);
      if (reaction->unit_i == 0) {
	fprintf(rxn_echo_fp,"DGZERO-UNITS\tkCal/mol\n");
      } else {
	fprintf(rxn_echo_fp,"DGZERO-UNITS\tkJ/mol\n");
      }
      fprintf(rxn_echo_fp,"//\n");
      reaction += 1; /* Caution address arithmetic. */
    }
    fclose(rxn_echo_fp);
  }
  return(success);
}
