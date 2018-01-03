/* size_file.c
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

#include "size_file.h"
int size_file(const char *path, int64_t *file_size) {
  /*
    Given the pathname of a file, determine its size in bytes using
    the stat function. return 1 on successfully setting *file_size else
    return 0.
    Called by: boltzmann_load.
    Calls:     stat, fprintf,strerror,fflush.

  */
  struct stat file_stat_buff;
  int ierr;
  int stat_err;
  int success;
  int padi;
  success = 1;
  ierr = stat(path,&file_stat_buff);
  stat_err  = errno; 
  if (ierr < 0) {
    fprintf(stderr,"size_file: Error, file %s not found/accesible.\n"
	    "Error was %d:%s\n",
	    path,stat_err,strerror(stat_err));
    fflush(stderr);
    *file_size = (int64_t)0;
    success = 0;
  } else {
    *file_size = (int64_t)file_stat_buff.st_size;
  }
  return(success);
}
