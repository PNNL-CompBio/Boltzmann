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

#include "boltzmann_structs.h"

#include "print_reactions_view.h"
int print_reactions_view(struct state_struct *state) {
  /*
    Print the reaction likelihoods per reaction to a file along with
    the reaction title and stoichiometric statement.
    Called by: boltzmann_run
    Calls:     fopen, fprintf, fclose (intrinsic)
  */
  struct rxn_struct *reactions;
  struct rxn_struct *reaction;
  struct rxn_matrix_struct *rxns_matrix;
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
  char *rxn_title_text;
  char *title;
  char *molecule;

  int success;
  int rxns;

  int ns;
  int padding;

  int np;
  int nr;
  
  int coeff;
  int j;

  int k;
  int rxn_view_hist_length;

  int nrxns;
  int net_fire;

  FILE *rxn_view_fp;
  FILE *lfp;

  char tab;
  char padc[7];

  tab  = (char)9;
  success = 1;
  molecules_text   = state->molecules_text;
  rxn_title_text   = state->rxn_title_text;

  rxn_view_fp = fopen(state->rxn_view_file,"w+");
  if (rxn_view_fp == NULL) {
    fprintf(stderr,
	    "print_reactions_view: Error could not open %s file.\n",
	    state->rxn_view_file);
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
    rxn_view_data  	 = state->rxn_view_likelihoods;
    rev_rxn_view_data 	 = state->rev_rxn_view_likelihoods;
    rxn_fire          	 = state->rxn_fire;
    rxn_view_hist_length = state->rxn_view_hist_length;
    nrxns                = (int)state->number_reactions;
    no_op_likelihood     = state->no_op_likelihood;
    /*
      Print out the number of times no reaction was chosen and its likelihoods.
      Skipping the net reaction fire and activities fields.
    */
    fprintf(rxn_view_fp," No Reaction\t0<=>0\t%d\t\t",rxn_fire[nrxns+nrxns]);
    for(k=0;k<rxn_view_hist_length;k++) {
      fprintf(rxn_view_fp,"\t1.0");
    }
    fprintf(rxn_view_fp,"\n");
    for (rxns=0;rxns < nrxns;rxns++) {
      if (reaction->title>=0) {
	title = (char *)&rxn_title_text[reaction->title];
	fprintf(rxn_view_fp,"%s\t",title);
      }
      nr = 0;
      for (j=rxn_ptrs[rxns];j<rxn_ptrs[rxns+1];j++) {
	coeff = coefficients[j];
	if (coeff < 0) {
	  if (nr > 0) {
	    /*
	      If this is not the first reactant printed out
	      put a + after the last one and before this one.
	    */
	    fprintf(rxn_view_fp," + ");
	  }
	  if (coeff < -1) {
	    coeff = -coeff;
	    fprintf(rxn_view_fp,"%d ",coeff);
	  }
	  molecule = (char*)&molecules_text[matrix_text[j]];
	  fprintf(rxn_view_fp,"%s",molecule);
	  nr += 1;
	}
      }
      fprintf(rxn_view_fp," => ");
      np = 0;
      for (j=rxn_ptrs[rxns];j<rxn_ptrs[rxns+1];j++) {
	coeff = coefficients[j];
	if (coeff > 0) {
	  if (np > 0) {
	    /*
	      If this is not the first product  printed out
	      put a + after the last one and before this one.
	    */
	    fprintf(rxn_view_fp," + ");
	  }
	  if (coeff > 1) {
	    fprintf(rxn_view_fp,"%d ",coeff);
	  }
	  molecule = (char*)&molecules_text[matrix_text[j]];
	  fprintf(rxn_view_fp,"%s",molecule);
	  np += 1;
	}
      }
      /*
	Print out the rxn fire count.
      */
      net_fire = rxn_fire[rxns] - rxn_fire[rxns+nrxns];
      fprintf(rxn_view_fp,"\t%d\t%d",rxn_fire[rxns],net_fire);
      /*
	Print out the rxn activity level.
      */
      fprintf(rxn_view_fp,"\t%le",activities[rxns]);
      /*
	Now print the reaction likelihoods for this reaction.
      */
      for(k=0;k<rxn_view_hist_length;k++) {
	fprintf(rxn_view_fp,"\t%le",*rxn_view_data);
	rxn_view_data++; /* Caution address arithmetic here.*/
      }
      fprintf(rxn_view_fp,"\n");
      /*
	Now repeat the process for the reverse reaction.
      */
      if (reaction->title >= 0) {
	title = (char *)&rxn_title_text[reaction->title];
	fprintf(rxn_view_fp,"%s\t",title);
      }
      nr = 0;
      for (j=rxn_ptrs[rxns];j<rxn_ptrs[rxns+1];j++) {
	coeff = coefficients[j];
	if (coeff < 0) {
	  if (nr > 1) {
	    fprintf(rxn_view_fp," + ");
	  }
	  if (coeff < -1) {
	    coeff = -coeff;
	    fprintf(rxn_view_fp,"%d ",coeff);
	  }
	  molecule = (char*)&molecules_text[matrix_text[j]];
	  fprintf(rxn_view_fp,"%s",molecule);
	  nr += 1;
	}
      }
      fprintf(rxn_view_fp," <= ");
      np = 0;
      for (j=rxn_ptrs[rxns];j<rxn_ptrs[rxns+1];j++) {
	coeff = coefficients[j];
	if (coeff > 0) {
	  if (np > 0) {
	    fprintf(rxn_view_fp," + ");
	  }
	  if (coeff > 1) {
	    fprintf(rxn_view_fp,"%d ",coeff);
	  }
	  molecule = (char*)&molecules_text[matrix_text[j]];
	  fprintf(rxn_view_fp,"%s",molecule);
	  np += 1;
	}
      }
      /*
	Print out the reverse rxn_fire count.
      */
      net_fire = -net_fire;
      fprintf(rxn_view_fp,"\t%d\t%d",rxn_fire[nrxns+rxns],net_fire);
      /*
	Print out the rxn activity level.
      */
      fprintf(rxn_view_fp,"\t%le",activities[rxns]);
      /*
	Now print the reaction likelihoods for this reverse reaction.
      */
      for(k=0;k<rxn_view_hist_length;k++) {
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
