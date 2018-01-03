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
    Allocate space for the species_name field for reading
    in intial concentrations, the concentrations vector,
    and the species matrix and its fields. The species 
    matrix is a transpose of the reactions matris.
    Called by: boltzmann
  */
  struct species_matrix_struct sms;
  struct species_matrix_struct *species_matrix;
  int64_t ask_for;
  int64_t one_l;
  int64_t usage;
  int64_t align_len;
  int64_t align_mask;
  int nu_species;
  int max_species_len;
  int success;
  int padi;
  int nzr;
  success = 1;
  one_l      = (int64_t)1;
  usage      = state->usage;
  align_mask = state->align_mask;
  align_len  = state->align_len;
  nu_species = state->unique_species;
  nzr        = state->number_species;
  max_species_len = state->max_species_len + 1;
  /*
    Allocate space for species name when reading initial 
    concentrations file.
  */
  ask_for    = max_species_len + 
    ((align_len - (max_species_len & align_mask)) & align_mask);
  usage += ask_for;
  state->species_name = (char *)calloc(one_l,ask_for);
  if (state->species_name == NULL) {
    fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	    "species_name field\n",ask_for);
    fflush(stderr);
    success = 0;
  }
  if (success) {
    /*
      Allocate space for the concentrations buffer.
    */
    ask_for = ((int64_t)nu_species) * ((int64_t)sizeof(double));
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
      Allocate space for the species_matrix;
    */
    ask_for = (int64_t)sizeof(sms);
    usage += ask_for;
    species_matrix = (struct species_matrix_struct *)calloc(one_l,ask_for);
    if (species_matrix == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "species_matrix field\n",ask_for);
      fflush(stderr);
      success = 0;
    } else {
      state->species_matrix = species_matrix;
    }
  }
  /*
    Allocate space for the fields of the species matrix.
  */
  if (success) {
    ask_for = ((int64_t)(nu_species+1)) * ((int64_t)sizeof(int64_t));
    usage += ask_for;
    species_matrix->species_ptrs = (int64_t*)calloc(one_l,ask_for);
    if (species_matrix->species_ptrs == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "species_ptrs field in species_matrix\n",ask_for);
      fflush(stderr);
      success = 0;
    } 
  }
  if (success) {
    ask_for = ((int64_t)(nzr)) * ((int64_t)sizeof(int64_t));
    usage += ask_for;
    species_matrix->rxn_indices = (int64_t*)calloc(one_l,ask_for);
    if (species_matrix->rxn_indices == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "rxn_indices field in species_matrix\n",ask_for);
      fflush(stderr);
      success = 0;
    } 
  }
  if (success) {
    ask_for = ((int64_t)(nzr)) * ((int64_t)sizeof(int64_t));
    usage += ask_for;
    species_matrix->rxn_indices = (int64_t*)calloc(one_l,ask_for);
    if (species_matrix->rxn_indices == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "rxn_indices field in species_matrix\n",ask_for);
      fflush(stderr);
      success = 0;
    } 
  }
  if (success) {
    ask_for = ((int64_t)(nzr)) * ((int64_t)sizeof(int64_t));
    usage += ask_for;
    species_matrix->coefficients = (int64_t*)calloc(one_l,ask_for);
    if (species_matrix->coefficients == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "coefficients field in species_matrix\n",ask_for);
      fflush(stderr);
      success = 0;
    } 
  }
  if (success) {
    ask_for = ((int64_t)(nu_species+1)) * sizeof(int64_t);
    state->transpose_work = (int64_t*)calloc(one_l,ask_for);
    if (species_matrix->coefficients == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "state->transpose_work field.\n",ask_for);
      fflush(stderr);
      success = 0;
    } 
  }
  state->usage = usage;
  return(success);
}
