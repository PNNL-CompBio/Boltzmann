/* boltzmann_load.c
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

#include "djb_timing_b.h"
#include "boltzmann_structs.h"
#include "size_file.h"
#include "boltzmann_mmap_superstate.h"
/*
#define DBG_BOLTZMANN_LOAD
*/
#include "boltzmann_load.h"
int boltzmann_load(char *superstate_filename, 
		   struct super_state_struct **super_statep) {
  /*
    Create a pointer to the superstate of a memorymapped superstate file,
    previously created by a call to boltzmann_boot.
    Called by: User (part of API)
    Calls:     size_file, boltzmann_mmap_superstate
  */
  int64_t superstate_size;
  int success;
  int padi;
  success = size_file(superstate_filename, &superstate_size);
  if (success) {
    success = boltzmann_mmap_superstate(superstate_filename,
					superstate_size,
					super_statep);
  } else {
    fprintf(stderr,"boltzmann_load unable to stat %s\n",superstate_filename);
    fflush(stderr);
    success = 0;
  }
  return(success);
}
