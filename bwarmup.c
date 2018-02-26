/* bwarmup.c
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

/*
  Maybe only for linux systems.
#ifdef LIBUNWIND
#include <libunwind.h>
#include <unwind.h>
#include "luwtb.h"
#endif
*/

#ifdef TIMING_ON
struct timing_struct timing_data;
#endif

/*
*/
#define BOLTZMANN_DBG 1
#include "boltzmann_init.h"
#include "bwarmup_run.h"
int main(int argc, char **argv)
{
  /*
    Calls:
      boltzmann_init
      boltzmann_run
  */
  struct state_struct *state;
  char *param_file_name;
  int success;
  int padi;
  /*
#ifdef LIBUNWIND
#include "luwtb1.h"
#endif
#ifdef LIBUNWIND
#include "luwtb2.h"
#endif
  */
  TIMING_INIT("./timingi.h");
  TIMING_START(TOTAL_TIME);
  TIMING_START(INITIALIZE);
  if (argc > 1) {
    param_file_name = argv[1];
  } else {
    param_file_name = NULL;
  }
  success = boltzmann_init(param_file_name,&state);
  TIMING_STOP(INITIALIZE);
  if (success) {
    state->print_ode_concs = 1;
    success = bwarmup_run(state);
  }
  TIMING_STOP(TOTAL_TIME);
  TIMING_PRINT(stdout);
  fflush(stdout);
  exit(0);
}
