/* boltzmann_load_aux_data.c
*******************************************************************************
Boltzmann

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
#include "boltzmann_flatten_aux_data.h"
#include "boltzmann_load_aux_data.h"
int boltzmann_load_aux_data(void **aux_data_p, struct state_struct *state) {
  /*
    set string pointers in the boltzman state struct 
    from memory pointed to by aux_data_p.
    Called by: biocellion/user, boltzman_test_save_load
    Calls: boltzmann_flatten_aux_data
  */
  int direction;
  int success;
  direction = 1;
  success = boltzmann_flatten_aux_data(state,aux_data_p,direction);
  return(success);
}
