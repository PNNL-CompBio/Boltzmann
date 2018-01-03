/* boltzmann_mmap_superstate.c
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

#include "boltzmann_mmap_superstate.h"

int boltzmann_mmap_superstate(char *global_state_filename,
			      int64_t mmap_file_len,
			      int64_t **super_statep) {
  /*
    mmap the super state file.
    Called by: boltzmann_boot, and user interface.
    Calls:     fopen, fileno, mmap, fprintf, fflush, strerror
  */
  int64_t *super_statev;
  int64_t mmap_start_addr;
  int64_t mmap_offset;
  int64_t zero_l;

  int mmap_err;
  
  int mmap_sm_flags;
  int mmap_read_prot;

  int success;
  int global_state_fd;
  FILE *global_state_fp;
  success         = 1;
  global_state_fp = fopen(global_state_filename,"r");
  if (global_state_fp == NULL) {
    fprintf(stderr,"boltzmann_mmap_superstate: Error, could not open "
	    "global state file for reading.%s\n",global_state_filename);
    fflush(stderr);
    success = 0;
  }
  if (success) {
    global_state_fd = fileno(global_state_fp);
    zero_l    = (int64_t)0;
    mmap_start_addr = zero_l;
    mmap_read_prot  = PROT_READ;
    mmap_sm_flags   = MAP_SHARED;
    mmap_offset     = zero_l;
    super_statev = (int64_t *) mmap(&mmap_start_addr,
				    mmap_file_len,
				    mmap_read_prot,
				    mmap_sm_flags,
				    global_state_fd,
				    mmap_offset);
    mmap_err = errno;
    if ((int64_t)(super_statev) < zero_l) {
      success = 0;
      fprintf(stderr,"boltzmann_mmap_superstate: Error unable to mmap "
	      "global state file, length = %ld, errno = %d:%s\n",
		mmap_file_len,mmap_err,strerror(mmap_err));
      fflush(stderr);
    } else {
      *super_statep = super_statev;
    }
  }
  return(success);
}
