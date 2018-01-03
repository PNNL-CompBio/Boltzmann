/* boltzmann_number_of_reaction_files.c
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

#include "boltzmann_number_of_reaction_files.h"
int boltzmann_number_of_reaction_files(struct super_state_struct *super_state) {
  /*
    Return the number of reaction files represented in the superstate.
    Called by : User (part of API).
  */
  int result;
  int padi;
  if (super_state) {
    result = (int)super_state->number_of_reaction_files;
  } else {
    result = 0;
  }
  return (result);
}
