/* print_reactions_view.c
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

#include "print_reactions_view.h"
int print_reactions_view(struct state_struct *state) {
  /*
    Print the reaction likelihoods per reaction to a file along with
    the reaction title and stoichiometric statement.
    Called by: boltzmann_run_sim
    Calls:     fopen, fprintf, fclose (intrinsic)
  */
  struct rxn_struct *reactions;
  struct rxn_struct *reaction;
  struct rxn_matrix_struct *rxns_matrix;
  int64_t *rxn_ptrs;
  int64_t *molecules_indices;
  int64_t *coefficients;
  char **matrix_text;
  double *rxn_view_data;
  double *rev_rxn_view_data;
  int    *rxn_fire;

  int success;
  int rxns;

  int ns;
  int padding;

  int np;
  int nr;
  
  int coeff;
  int j;

  int k;
  int lthl;

  int nrxns;
  int net_fire;

  FILE *rxn_view_fp;
  FILE *lfp;

  char tab;
  char padc[7];

  tab  = (char)9;
  success = 1;

  rxn_view_fp = fopen(state->rxn_view_file,"w+");
  if (rxn_view_fp == NULL) {
    fprintf(stderr,
	    "print_reactions_view: Error could not open %s file.\n",
	    state->rxn_view_file);
    success = 0;
  }
  if (success) {
    reaction       = state->reactions;
    rxns_matrix    = state->reactions_matrix;
    rxn_ptrs       = rxns_matrix->rxn_ptrs;
    molecules_indices = rxns_matrix->molecules_indices;
    coefficients   = rxns_matrix->coefficients;
    matrix_text    = rxns_matrix->text;
    rxn_view_data  = state->rxn_view_likelihoods;
    rev_rxn_view_data = state->rev_rxn_view_likelihoods;
    rxn_fire          = state->rxn_fire;
    lthl              = state->lthl;
    nrxns             = state->number_reactions;
    for (rxns=0;rxns < nrxns;rxns++) {
      if (reaction->title) {
	fprintf(rxn_view_fp,"%s\t",reaction->title);
      }
      nr = 0;
      for (j=rxn_ptrs[rxns];j<rxn_ptrs[rxns+1];j++) {
	coeff = coefficients[j];
	if (coeff < 0) {
	  if (coeff < -1) {
	    coeff = -coeff;
	    fprintf(rxn_view_fp,"%d ",coeff);
	  }
	  fprintf(rxn_view_fp,"%s",matrix_text[j]);
	  nr += 1;
	  if (nr < reaction->num_reactants) {
	    fprintf(rxn_view_fp," + ");
	  } else {
	    fprintf(rxn_view_fp," => ");
	    break;
	  }
	}
      }
      np = 0;
      for (j=rxn_ptrs[rxns];j<rxn_ptrs[rxns+1];j++) {
	coeff = coefficients[j];
	if (coeff > 0) {
	  if (coeff > 1) {
	    fprintf(rxn_view_fp,"%d ",coeff);
	  }
	  fprintf(rxn_view_fp,"%s",matrix_text[j]);
	  np += 1;
	  if (np < reaction->num_products) {
	    fprintf(rxn_view_fp," + ");
	  } else {
	    fprintf(rxn_view_fp," ");
	    break;
	  }
	}
      }
      /*
	Print out the rxn fire count.
      */
      net_fire = rxn_fire[rxns] - rxn_fire[rxns+nrxns];
      fprintf(rxn_view_fp,"\t%d\t%d",rxn_fire[rxns],net_fire);
      /*
	Now print the reaction likelihoods for this reaction.
      */
      for(k=0;k<lthl;k++) {
	fprintf(rxn_view_fp,"\t%le",*rxn_view_data);
	rxn_view_data++; /* Caution address arithmetic here.*/
      }
      fprintf(rxn_view_fp,"\n");
      /*
	Now repeat the process for the reverse reaction.
      */
      if (reaction->title) {
	fprintf(rxn_view_fp,"%s\t",reaction->title);
      }
      nr = 0;
      for (j=rxn_ptrs[rxns];j<rxn_ptrs[rxns+1];j++) {
	coeff = coefficients[j];
	if (coeff < 0) {
	  if (coeff < -1) {
	    coeff = -coeff;
	    fprintf(rxn_view_fp,"%d ",coeff);
	  }
	  fprintf(rxn_view_fp,"%s",matrix_text[j]);
	  nr += 1;
	  if (nr < reaction->num_reactants) {
	    fprintf(rxn_view_fp," + ");
	  } else {
	    fprintf(rxn_view_fp," <= ");
	    break;
	  }
	}
      }
      np = 0;
      for (j=rxn_ptrs[rxns];j<rxn_ptrs[rxns+1];j++) {
	coeff = coefficients[j];
	if (coeff > 0) {
	  if (coeff > 1) {
	    fprintf(rxn_view_fp,"%d ",coeff);
	  }
	  fprintf(rxn_view_fp,"%s",matrix_text[j]);
	  np += 1;
	  if (np < reaction->num_products) {
	    fprintf(rxn_view_fp," + ");
	  } else {
	    fprintf(rxn_view_fp," ");
	    break;
	  }
	}
      }
      /*
	Print out the reverse rxn_fire count.
      */
      net_fire = -net_fire;
      fprintf(rxn_view_fp,"\t%d\t%d",rxn_fire[nrxns+rxns],net_fire);
      /*
	Now print the reaction likelihoods for this reverse reaction.
      */
      for(k=0;k<lthl;k++) {
	fprintf(rxn_view_fp,"\t%le",*rev_rxn_view_data);
	rev_rxn_view_data++; /* Caution address arithmetic here.*/
      }
      fprintf(rxn_view_fp,"\n");
      reaction += 1; /* Caution address arithmetic here.*/
    }
    fclose(rxn_view_fp);
  }
  return(success);
}
