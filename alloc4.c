/* alloc4.c
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

#include "alloc4.h"

int alloc4(struct state_struct *state,
	   struct molecules_matrix_struct **molecules_matrix_p,
	   int64_t **transpose_work_p) {
  /*
    Allocate space for the molecules_matrix structure and its
    subfields and for the transpose_work space needed to form it.
    Called by: rxn_map_init
    Calls:     calloc, fprintf, fflush,
  */

  struct molecules_matrix_struct *molecules_matrix;
  struct molecules_matrix_struct sms;
  int64_t *transpose_work;
  int64_t usage;
  int64_t ask_for;
  int64_t one_l;
  int success;
  int nu_molecules;
  int nzr;
  one_l         = (int64_t)1;
  usage         = state->usage;
  nu_molecules  = (int)state->nunique_molecules;
  nzr           = (int)state->number_molecules;
  success       = 1;
  /*  
    Allocate space for the molecules_matrix;
  */
  if (success) {
    ask_for = (int64_t)sizeof(sms);
    usage += ask_for;
    molecules_matrix = (struct molecules_matrix_struct *)calloc(one_l,ask_for);
    if (molecules_matrix == NULL) {
      fprintf(stderr,"alloc4: Error unable to allocate %ld bytes for "
	      "molecules_matrix field\n",ask_for);
      fflush(stderr);
      success = 0;
    } else {
      *molecules_matrix_p = molecules_matrix;
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
      fprintf(stderr,"alloc4: Error unable to allocate %ld bytes for "
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
      fprintf(stderr,"alloc4: Error unable to allocate %ld bytes for "
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
      fprintf(stderr,"alloc4: Error unable to allocate %ld bytes for "
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
      fprintf(stderr,"alloc4: Error unable to allocate %ld bytes for "
	      "coefficients field in molecules_matrix\n",ask_for);
      fflush(stderr);
      success = 0;
    } 
  }
  if (success) {
    ask_for = ((int64_t)(nu_molecules+1)) * ((int64_t)sizeof(int64_t));
    usage += ask_for;
    transpose_work = (int64_t*)calloc(one_l,ask_for);
    if (molecules_matrix->coefficients == NULL) {
      fprintf(stderr,"alloc4: Error unable to allocate %ld bytes for "
	      "transpose_work vector.\n",ask_for);
      fflush(stderr);
      success = 0;
    } else {
      *transpose_work_p = transpose_work;
    }
  }
  state->usage = usage;
  return(success);
}
