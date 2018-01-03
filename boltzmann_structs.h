/* boltzmann_structs.h 
  
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
*/
#ifndef __BOLTZMANN_STRUCTS__
#define __BOLTZMANN_STRUCTS__ 1
/*
  System include files.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <float.h>
#include <signal.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
/*
  Maybe only for linux systems.
*/
#ifdef LIBUNWIND
#include <libunwind.h>
#include <unwind.h>
#endif

#include "djb_timing_b.h"
/*
  Data structures include files.
*/
#include "super_state_struct.h"
#include "super_state_pointers_struct.h"
#include "state_struct.h"
#include "rxn_struct.h"
#include "rxn_matrix_struct.h"
#include "molecules_matrix_struct.h"
#include "istring_elem_struct.h"
#include "vgrng_state_struct.h"
#endif
