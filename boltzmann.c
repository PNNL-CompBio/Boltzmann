/* boltzmann.c
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
*/
#ifdef LIBUNWIND
#include <libunwind.h>
#include <unwind.h>
#endif

#ifdef TIMING_ON
struct timing_struct timing_data;
#endif

#ifdef LIBUNWIND
#include "luwtb.h"
#endif
/*
*/
#define BOLTZMANN_DBG 1
#include "boltzmann_init.h"
#include "boltzmann_build_agent_data_block.h"
#include "boltzmann_run.h"
#include "print_boundary_flux.h"
#include "print_restart_file.h"
int main(int argc, char **argv)
{
  /*
    Calls:
      boltzmann_init
      boltzmann_run
  */
  struct state_struct *state;
  void *agent_data;
  char *param_file_name;
  int success;
  int padi;

#ifdef LIBUNWIND
#include "luwtb1.h"
#endif
#ifdef LIBUNWIND
#include "luwtb2.h"
#endif
  TIMING_INIT("./timingi.h");
  TIMING_START(TOTAL_TIME);
  TIMING_START(INITIALIZE);
  if (argc > 1) {
    param_file_name = argv[1];
  } else {
    param_file_name = NULL;
  }
  success = boltzmann_init(param_file_name,&state);
  if (success) {
    boltzmann_build_agent_data_block(state,&agent_data);
  }
  TIMING_STOP(INITIALIZE);
  if (success) {
    success = boltzmann_run(state,agent_data);
  }
  if (state->print_output == 0) {
    if (state->num_fixed_concs > (int64_t)0) {
      state->bndry_flux_fp = fopen(state->bndry_flux_file,"w");
      print_boundary_flux(state);
    } /* end if (state->num_fixed_concs ...) */
    if (success) {
      success = print_restart_file(state);
    }
  }
  TIMING_STOP(TOTAL_TIME);
  TIMING_PRINT(stdout);
  fflush(stdout);
  exit(0);
}
