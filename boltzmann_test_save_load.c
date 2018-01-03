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
#define BOLTZMANN_DBG 1
*/
#include "boltzmann_init.h"
#include "boltzmann_save_state.h"
#include "boltzmann_save_aux_data.h"
#include "boltzmann_load_state.h"
#include "boltzmann_load_aux_data.h"
#include "open_output_files.h"
#include "boltzmann_run.h"
int main(int argc, char **argv) {
  /*
    Calls:
      boltzmann_init
      boltzmann_run
  */
  struct state_struct *state;
  struct state_struct *recovered_state;
  void *flattened_state;
  char *param_file_name;
  char *aux_data;
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
    flattened_state = NULL;
    success = boltzmann_save_state(state,&flattened_state);
  }
  if (success) {
    success = boltzmann_save_aux_data(state,(void**)&aux_data);
  }
  if (success) {
    boltzmann_load_state(flattened_state,&recovered_state);
    recovered_state->lfp = state->lfp;
  }
  if (success) {
    success = boltzmann_load_aux_data((void**)&aux_data,recovered_state);
  }
  if (success) {
    if (recovered_state->print_output) {
      /*
	If you were not just reloading a just saved state and you wanted
        printing turned on, you would reopen the output files.
      open_output_files(recovered_state);
      As it is we will just copy file pointers.
      */
      recovered_state->rxn_fp         = state->rxn_fp;
      recovered_state->conc_fp        = state->conc_fp;
      recovered_state->out_fp         = state->out_fp;
      recovered_state->counts_out_fp  = state->counts_out_fp;
      recovered_state->concs_out_fp   = state->concs_out_fp;
      recovered_state->ode_concs_fp   = state->ode_concs_fp;
      recovered_state->net_lklhd_fp   = state->net_lklhd_fp;
      recovered_state->free_energy_fp = state->free_energy_fp;
      recovered_state->restart_fp     = state->restart_fp;
      recovered_state->rxn_view_fp    = state->rxn_view_fp;
      recovered_state->bndry_flux_fp  = state->bndry_flux_fp;
      recovered_state->cmpt_fp        = state->cmpt_fp;
      recovered_state->ode_dconcs_fp  = state->ode_dconcs_fp;
      recovered_state->ode_lklhd_fp   = state->ode_lklhd_fp;
      recovered_state->ode_bflux_fp   = state->ode_bflux_fp;
    }
  }
  TIMING_STOP(INITIALIZE);
  if (success) {
    success = boltzmann_run(recovered_state);
  }
  TIMING_STOP(TOTAL_TIME);
  TIMING_PRINT(stdout);
  fflush(stdout);
  exit(0);
}
