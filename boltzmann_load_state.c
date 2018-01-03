/* boltzmann_load_state.c
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
#include "boltzmann_flatten_state.h"
#include "boltzmann_build_agent_data_block.h"
#include "boltzmann_load_state.h"
int boltzmann_load_state(void *flattened_state,
			 struct state_struct **state_p, 
			 void **agent_data_p) {
  /*
    Allocate and load a boltzmann state from a flattened_state
    created by boltzman_save_state, also allocate and fill
    an agent_data_block 
    Called by: biocellion/user, boltzmann_test_save_load
    Calls:     boltzmann_flatten_state
  */
  struct state_struct *state;
  void *in_flattened_state;
  int direction;
  int success;
  in_flattened_state = flattened_state;
  direction = 1;
  state = NULL;
  success = boltzmann_flatten_state(&state,&in_flattened_state,direction,
				    stderr);
  if (success) {
    success = boltzmann_build_agent_data_block(state,agent_data_p);
  }
  *state_p = state;
  return(success);
}
