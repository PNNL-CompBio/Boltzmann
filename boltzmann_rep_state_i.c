/* boltzmann_rep_state_i.c
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


#include "djb_timing_b.h"
#include "boltzmann_structs.h"

#include "flatten_state.h"
#include "boltzmann_rep_state_i.h"
int boltzmann_rep_state_i(int64_t *super_statep, int i, 
			 struct state_struct **l_statep) {
  struct  state_struct *lstate;
  int64_t *super_statel;
  int64_t *state_offsets_sizes;
  char    *super_statec;
  int64_t lstate_offset;
  int64_t lstate_size;
  int64_t one_l;
  int success;
  int index;
  success = 1;
  one_l   = (int64_t)1;
  super_statec = (char*)super_statep;
  super_statel = (int64_t*)super_statep;
  state_offsets_sizes = &super_statel[super_statel[12]];
  if ((int64_t)i < super_statel[0]) {
    index = i<<1;
    lstate_offset = state_offsets_sizes[index];
    lstate_size   = state_offsets_sizes[index+1];
    lstate        = (struct state_struct*)calloc(one_l,lstate_size);
    if (lstate) {
      memcpy(lstate,(void *)&super_statec[lstate_offset],lstate_size);
      success = flatten_state(lstate,&lstate);
    } else {
      success = 0;
    }
  } else {
    lstate = NULL;
    success = 0;
  }
  *l_statep = lstate;
  return(success);
}
