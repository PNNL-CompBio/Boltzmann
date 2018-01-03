/* boltzmann_length_state_i.c
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
#include <errno.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "boltzmann_structs.h"

#include "boltzmann_length_state_i.h"
int64_t boltzmann_length_state_i(struct super_state_struct *super_state,
				 int state_index) {
  /*
    Return the length of local state i, if it exists.
    Called by User (part of the API).
  */
  int64_t *super_statel;
  int64_t *state_offsets_sizes;
  int64_t number_of_reaction_files;
  int64_t state_offsets_sizes_offset;
  int64_t length;
  if (super_state) {
    number_of_reaction_files = super_state->number_of_reaction_files;
    if (state_index < number_of_reaction_files) {
      state_offsets_sizes_offset = super_state->state_offsets_sizes_offset_in_words;
      super_statel = (int64_t *)super_state;
      state_offsets_sizes = (int64_t *)&super_statel[state_offsets_sizes_offset];
      length = state_offsets_sizes[(state_index<<1)+1];
    } else {
      length = (int64_t)0;
    }
  } else {
    length = (int64_t)0;
  }
  return(length);
}
