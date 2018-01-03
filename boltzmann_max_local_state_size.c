/* boltzmann_max_local_state_size.c
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
/*
#define DBG_BOLTZMANN_MAX_LOCAL_STATE_SIZE 1
*/
#include "boltzmann_max_local_state_size.h"
int64_t boltzmann_max_local_state_size(struct super_state_struct *super_state) {
  int64_t size;
  if (super_state) {
    size = super_state->maximum_state_size_in_bytes;
  } else {
    size = (int64_t)0;
  }
  return(size);
}
