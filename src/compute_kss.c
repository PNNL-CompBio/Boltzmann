/* compute_kss.c
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
 *  Created on: Apr 21, 2014
 *      Author:  Doug Baxter
*/
#include "boltzmann_structs.h"
#include "conc_to_pow.h"
#include "compute_kss.h"

int compute_kss(struct state_struct *state) {
  /*
    Set the delta_g0 field for all of the reactions
    based on energies of formation.
    Called by: energy_init
    Calls:     conc_to_pow, fprintf, fflush
  */
  struct reaction_struct *reactions;
  struct reaction_struct *reaction;

  struct reactions_matrix_struct *rxns_matrix;
  struct compartment_struct *compartments;
  struct compartment_struct *compartment;
  double *kss_e_val;
  double *kss_u_val;

  double *kss;
  double *kssr;
  double *ke;
  double *coefficients;

  double coeff;
  double factorial;

  int64_t *rxn_ptrs;
  int64_t *molecules_indices;
  int64_t *compartment_indices;
  int64_t *matrix_text;

  double  rrxn_kss;
  double  kss_r_mod;
  double  kss_p_mod;
  double  ntotal_exp;
  double  ntotal_opt;
  double  conc_to_count;
  double  count_to_conc;
  int64_t print_output;

  char   *molecules_text;

  int nrxns;
  int rxns;

  int j;
  int k;

  int ci;
  int success;

  FILE* lfp;
  FILE* efp;

  success = 1;
  
  reactions       	   = state->reactions;
  reaction                 = reactions;
  rxns_matrix    	   = state->reactions_matrix;
  print_output             = state->print_output;
  rxn_ptrs       	   = rxns_matrix->rxn_ptrs;
  molecules_indices        = rxns_matrix->molecules_indices;
  compartment_indices      = rxns_matrix->compartment_indices;
  coefficients   	   = rxns_matrix->coefficients;
  matrix_text    	   = rxns_matrix->text;
  nrxns                    = (int)state->number_reactions;
  molecules_text           = state->molecules_text;
  kss                      = state->kss;
  kssr                     = state->kssr;
  kss_e_val                = state->kss_e_val;
  kss_u_val                = state->kss_u_val;
  ke                       = state->ke;
  compartments             = state->sorted_compartments;
  lfp                      = state->lfp;
  print_output             = print_output && lfp;
  if (print_output) {
    fprintf(lfp, "Output from compute_kss.c: \n");
    fprintf(lfp,"compute_kss: number of reactions = %i\n",nrxns);
  }
  if (state->adjust_steady_state) {
    for (rxns=0;rxns < nrxns;rxns++) {
      if (print_output) {
        fprintf(lfp,"compute_kss: reaction no. %i\n",rxns);
        fprintf(lfp,"compute_kss: number of reactants = %i\n",
    	      reaction->num_reactants);
        fprintf(lfp,"compute_kss: Number of products = %i\n",
    	      reaction->num_products);
      }
      /*
        Set the steady state adjustment value for the reaction.
      */
      rrxn_kss = 1.0;
      for (j=rxn_ptrs[rxns];j<rxn_ptrs[rxns+1];j++) {
        k = molecules_indices[j];
	ci = compartment_indices[j];
	compartment = (struct compartment_struct *)&compartments[ci];
	ntotal_exp               = compartment->ntotal_exp;
	ntotal_opt               = compartment->ntotal_opt;
        coeff = coefficients[j];
	conc_to_count = compartment->conc_to_count;
	count_to_conc = compartment->count_to_conc;
        if (coeff < 0) {
	  coeff = -coeff;
	  kss_r_mod = kss_e_val[k] * conc_to_count;
	  factorial = 0.0;
	  rrxn_kss = rrxn_kss * conc_to_pow(kss_r_mod,coeff,factorial);
	  /*
	  for (i=0;i<coeff;i++) {
	    rrxn_kss = rrxn_kss * kss_r_mod;
	  }
	  */
        } else {
	  if (coeff > 0) {
	    if (kss_e_val[k] > 0.0) {
	      kss_p_mod = count_to_conc/kss_e_val[k];
	      factorial = 0.0;
	      rrxn_kss = rrxn_kss * conc_to_pow(kss_p_mod,coeff,factorial);
	      /*
	      for (i=0;i<coeff;i++) {
		rrxn_kss = rrxn_kss * kss_p_mod;
	      }
	      */
	    } else {
	      /*
		Question here is should we set rxn_kss to 1, or 0,
		or just not modify it at all?
	      */
	    }
	  }
	}
      }
      rrxn_kss = ke[rxns] * rrxn_kss;
      kssr[rxns] = rrxn_kss;
      if (rrxn_kss > 0.0) {
	kss[rxns] = 1.0/rrxn_kss;
      } else {
	kss[rxns] = 1.0;
	kssr[rxns] = 1.0;
      }
      if (print_output) {
        fprintf(lfp,"Computed kss = %le\n",kss[rxns]);
      }
      reaction += 1; /* Caution address arithmetic */
    } 
  } else {
    for (rxns=0;rxns < nrxns;rxns++) {
      kss[rxns] = 1.0;
      kssr[rxns] = 1.0;
    }
  }
  if (print_output) {
    fprintf(lfp,"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
  }
  return (success);
}
