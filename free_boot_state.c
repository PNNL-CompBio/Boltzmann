/* free_boot_state.c
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

#include "djb_timing_b.h"
#include "boltzmann_structs.h"
/*
#define DBG_FREE_BOOT_STATE 1
*/
#include "free_boot_state.h"
int free_boot_state(struct state_struct **statep) {
  /*
    Called by: boltzmann_init
    Calls:     free.
  */
  struct state_struct *state;
  struct rxn_matrix_struct *reactions_matrix; 
  int success;
  int padi;
  success = 1;
  state = *statep;
  if (state) {
    if (state->params_file) {
      free(state->params_file);
    }
    if (state->param_buffer) {
      free(state->param_buffer);
    }
    if (state->rxn_file_keyword_buffer) {
      free(state->rxn_file_keyword_buffer);
    }
    if (state->rxn_file_keywords) {
      free(state->rxn_file_keywords);
    }
    if (state->rxn_file_keyword_lengths) {
      free(state->rxn_file_keyword_lengths);
    }
    if (state->rxn_title_text) {
      free(state->rxn_title_text);
    }
    if (state->reactions) {
      free(state->reactions);
    }
    reactions_matrix = state->reactions_matrix;
    if (reactions_matrix) {
      if (reactions_matrix->rxn_ptrs) {
	free(reactions_matrix->rxn_ptrs);
      }
      if (reactions_matrix->molecules_indices) {
	free(reactions_matrix->molecules_indices);
      }
      if (reactions_matrix->compartment_indices) {
	free(reactions_matrix->compartment_indices);
      }
      if (reactions_matrix->coefficients) {
	free(reactions_matrix->coefficients);
      }
      if (reactions_matrix->text) {
	free(reactions_matrix->text);
      }
      free(state->reactions_matrix);    
    }
    if (state->unsorted_molecules) {
      free(state->unsorted_molecules);
    }
    if (state->unsorted_cmpts) {
      free(state->unsorted_cmpts);
    }
    if (state->activities) {
      free(state->activities);
    }
    if (state->vgrng_state) {
      free(state->vgrng_state);
    }
    if (state->vgrng2_state) {
      free(state->vgrng2_state);
    }
    if (state->compartment_ptrs) {
      free(state->compartment_ptrs);
    }
    if (state->current_concentrations) {
      free(state->current_concentrations);
    }
    if (state->bndry_flux_concs) {
      free(state->bndry_flux_concs);
    }
    if (state->dg0s) {
      free(state->dg0s);
    }
    if (state->ke) {
      free(state->ke);
    }
    if (state->workspace_base) {
      free(state->workspace_base);
    }
    free(state);
  }
  return(success);
}

