/* alloc3.c
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

#include "boltzmann_structs.h"

#include "alloc3.h"

int alloc3(struct state_struct *state) {
  /*
    Allocate space for the molecule_name field for reading
    in intial concentrations, the concentrations vector,
    and the molecules matrix and its fields. The molecules 
    matrix is a transpose of the reactions matris.
    Called by: boltzmann
  */
  struct molecules_matrix_struct sms;
  struct molecules_matrix_struct *molecules_matrix;
  int64_t ask_for;
  int64_t one_l;
  int64_t usage;
  int64_t align_len;
  int64_t align_mask;
  int nu_molecules;
  int max_molecule_len;

  int success;
  int nrxns;

  int nzr;
  int padi;
  success = 1;
  one_l      = (int64_t)1;
  usage      = state->usage;
  align_mask = state->align_mask;
  align_len  = state->align_len;
  nu_molecules = state->unique_molecules;
  nzr        = state->number_molecules;
  nrxns      = state->number_reactions;
  max_molecule_len = state->max_molecule_len + 1;
  /*
    Allocate space for molecules name when reading initial 
    concentrations file.
  */
  ask_for    = max_molecule_len + 
    ((align_len - (max_molecule_len & align_mask)) & align_mask);
  usage += ask_for;
  state->molecule_name = (char *)calloc(one_l,ask_for);
  if (state->molecule_name == NULL) {
    fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	    "molecule_name field\n",ask_for);
    fflush(stderr);
    success = 0;
  }
  if (success) {
    /*
      Allocate space for the concentrations buffer.
    */
    ask_for = ((int64_t)nu_molecules) * ((int64_t)sizeof(double));
    usage += ask_for;
    state->concentrations = (double *)calloc(one_l,ask_for);
    if (state->concentrations == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "concentrations field\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  if (success) {
    /*
      Allocate space for the molecules_matrix;
    */
    ask_for = (int64_t)sizeof(sms);
    usage += ask_for;
    molecules_matrix = (struct molecules_matrix_struct *)calloc(one_l,ask_for);
    if (molecules_matrix == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "molecules_matrix field\n",ask_for);
      fflush(stderr);
      success = 0;
    } else {
      state->molecules_matrix = molecules_matrix;
    }
  }
  /*
    Allocate space for the fields of the molecules matrix.
  */
  if (success) {
    ask_for = ((int64_t)(nu_molecules+1)) * ((int64_t)sizeof(int64_t));
    usage += ask_for;
    molecules_matrix->molecules_ptrs = (int64_t*)calloc(one_l,ask_for);
    if (molecules_matrix->molecules_ptrs == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "molecules_ptrs field in molecules_matrix\n",ask_for);
      fflush(stderr);
      success = 0;
    } 
  }
  if (success) {
    ask_for = ((int64_t)(nzr)) * ((int64_t)sizeof(int64_t));
    usage += ask_for;
    molecules_matrix->rxn_indices = (int64_t*)calloc(one_l,ask_for);
    if (molecules_matrix->rxn_indices == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "rxn_indices field in molecules_matrix\n",ask_for);
      fflush(stderr);
      success = 0;
    } 
  }
  if (success) {
    ask_for = ((int64_t)(nzr)) * ((int64_t)sizeof(int64_t));
    usage += ask_for;
    molecules_matrix->rxn_indices = (int64_t*)calloc(one_l,ask_for);
    if (molecules_matrix->rxn_indices == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "rxn_indices field in molecules_matrix\n",ask_for);
      fflush(stderr);
      success = 0;
    } 
  }
  if (success) {
    ask_for = ((int64_t)(nzr)) * ((int64_t)sizeof(int64_t));
    usage += ask_for;
    molecules_matrix->coefficients = (int64_t*)calloc(one_l,ask_for);
    if (molecules_matrix->coefficients == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "coefficients field in molecules_matrix\n",ask_for);
      fflush(stderr);
      success = 0;
    } 
  }
  if (success) {
    ask_for = ((int64_t)(nu_molecules+1)) * ((int64_t)sizeof(int64_t));
    usage += ask_for;
    state->transpose_work = (int64_t*)calloc(one_l,ask_for);
    if (molecules_matrix->coefficients == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "state->transpose_work field.\n",ask_for);
      fflush(stderr);
      success = 0;
    } 
  }
  if (success) {
    ask_for = ((int64_t)nrxns) * ((int64_t)sizeof(double));
    usage += ask_for;
    state->dg0s = (double *)calloc(one_l,ask_for);
    if (state->dg0s == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "state->dg0s field.\n",ask_for);
      fflush(stderr);
      success = 0;
    } 
    
  }
  if (success) {
    ask_for = ((int64_t)nrxns) * ((int64_t)sizeof(double));
    usage += ask_for;
    state->dgs = (double *)calloc(one_l,ask_for);
    if (state->dgs == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "state->dgs field.\n",ask_for);
      fflush(stderr);
      success = 0;
    } 
    
  }
  if (success) {
    ask_for = ((int64_t)nrxns) * ((int64_t)sizeof(double));
    usage += ask_for;
    state->q_r = (double *)calloc(one_l,ask_for);
    if (state->q_r == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "state->q_r field.\n",ask_for);
      fflush(stderr);
      success = 0;
    } 
    
  }
  if (success) {
    ask_for = ((int64_t)nrxns) * ((int64_t)sizeof(double));
    usage += ask_for;
    state->ke = (double *)calloc(one_l,ask_for);
    if (state->q_r == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "state->ke field.\n",ask_for);
      fflush(stderr);
      success = 0;
    } 
  }
  state->usage = usage;
  return(success);
}
