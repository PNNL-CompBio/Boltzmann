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
#include "boltzmann_structs.h"
/*
#define DBG_FREE_BOOT_STATE 1
*/
#include "free_boot_state.h"
#include "free_boot_state2.h"
int free_boot_state(struct state_struct **statep) {
  /*
    Free temporary space used in setup.
    Called by: boltzmann_init
    Calls:     free_boot_state2, free.
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
    if (state->solvent_string) {
      free(state->solvent_string);
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
    if (state->workspace_base) {
      free(state->workspace_base);
    }
    success = free_boot_state2(statep);
    free(state);
  }
  return(success);
}

