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

#include "compute_kss.h"

int compute_kss(struct state_struct *state) {
  /*
    Set the delta_g0 field for all of the reactions
    based on energies of formation.
    Called by: boltzmann_init, boltzmann_boot
    Calls:     fprintf, fflush
  */
  int success;
  int64_t sum;


  double *kss_e_val;
  double *kss_u_val;

  double *kss;
  double *kssr;

  struct rxn_struct *reactions;
  struct rxn_struct *reaction;

  struct rxn_matrix_struct *rxns_matrix;
  struct molecule_struct *compartments;
  struct molecule_struct *compartment;
  int64_t *rxn_ptrs;
  int64_t *molecules_indices;
  int64_t *compartment_indices;
  int64_t *coefficients;
  int64_t *matrix_text;

  double  rxn_kss;
  double  kss_r_mod;
  double  kss_p_mod;
  double  ntotal_exp;
  double  ntotal_opt;
  int64_t print_output;
  int64_t mto;

  char   *molecule;
  char   *molecules_text;

  int nrxns;
  int rxns;

  int nr;
  int j;

  int k;
  int coeff;

  int i;
  int ci;

  FILE* lfp;


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
  compartments             = state->sorted_cmpts;
  lfp                  = state->lfp;
  print_output         = print_output && lfp;
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
      rxn_kss = 1.0;
      for (j=rxn_ptrs[rxns];j<rxn_ptrs[rxns+1];j++) {
        k = molecules_indices[j];
	ci = compartment_indices[j];
	compartment = (struct molecule_struct *)&compartments[ci];
	ntotal_exp               = compartment->ntotal_exp;
	ntotal_opt               = compartment->ntotal_opt;
        coeff = coefficients[j];
        if (coeff < 0) {
	  coeff = -coeff;
	  if (kss_u_val[k] > 0.0) {
	    kss_r_mod = (kss_e_val[k]/kss_u_val[k]) * (ntotal_opt/ntotal_exp);
	    for (i=0;i<coeff;i++) {
	      rxn_kss = rxn_kss * kss_r_mod;
	    }
	  } else {
	    /*
	      Question here is should we set rxn_kss to 1, or 0,
	      or just not modify it at all?
	    */
	  }
        } else {
	  if (coeff > 0) {
	    if (kss_e_val[k] > 0.0) {
	      kss_p_mod = (kss_u_val[k]/kss_e_val[k]) * (ntotal_exp/ntotal_exp);
	      for (i=0;i<coeff;i++) {
		rxn_kss = rxn_kss * kss_p_mod;
	      }
	    } else {
	      /*
		Question here is should we set rxn_kss to 1, or 0,
		or just not modify it at all?
	      */
	    }
	  }
	}
      }
      kss[rxns] = rxn_kss;
      if (rxn_kss > 0) {
	kssr[rxns] = 1.0/rxn_kss;
      } else {
	kssr[rxns] = 1.0;
      }
      if (print_output) {
        fprintf(lfp,"Computed kss = %le\n",rxn_kss);
      }
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
