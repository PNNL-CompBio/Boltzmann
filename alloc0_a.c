/* alloc0.c
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
#include "boltzmann_set_filename_ptrs.h"

#include "alloc0_a.h"
int alloc0_a(struct state_struct *state) {
  /*
    Allocate auxilliary data fields (fields use in setup and printing,
    but not needed for boltzmann_run) that are allocated in alloc0

    Called by: alloc0,
    Calls:     boltzmann_set_filename_ptrs
               fprintf,fflush
    Sets: The following fields in state:
       num_files,
       max_filename_len,

       params_file,
       reaction_file,
       init_conc_file,
       input_dir,
       output_file,
       log_file,
       counts_out_file,
       ode_concs_file,
       net_lklhd_file     
       nl_bndry_flx_file
       rxn_lklhd_file,
       free_energy_file,
       restart_file,
       rxn_view_file,
       bndry_flux_file,
       pseudoisomer_file,
       compartment_file,
       sbml_file,
       ms2js_file,
       kg2js_file,
       rxn_echo_file,
       rxn_mat_file,
       dg0ke_file,
       dictionary_file,
       ode_dconcs_file,
       ode_lklhd_file,
       ode_bflux_file,
       concs_out_file,
       ode_concs_file,

       solvent_string
  */
  int64_t max_file_name_len;
  int64_t num_state_files;
  int64_t one_l;
  int64_t ask_for;
  int64_t usage;
  char *solvent_string;
  int success;
  success = 1;
  max_file_name_len = (int64_t)128;
  num_state_files   = (int64_t)32;
  one_l             = (int64_t)1;
  usage             = state->usage;
  state->num_files        =  num_state_files;
  state->max_filename_len = max_file_name_len;
  ask_for = ((int64_t)num_state_files) * max_file_name_len;
  usage   += ask_for;
  state->params_file = (char *)calloc(one_l,ask_for);
  if (state->params_file == NULL) {
    success = 0;
    fprintf(stderr,
	    "alloc0_a: unable to allocate %ld bytes for state file names.\n",
	    ask_for);
    fflush(stderr);
  }
  if (success) {
    state->max_filename_len   = max_file_name_len;
    boltzmann_set_filename_ptrs(state);
    ask_for = ((int64_t)64) * sizeof(char);
    usage += ask_for;
    solvent_string = (char *)calloc(one_l,ask_for);
    if (solvent_string) {
      state->solvent_string = solvent_string;
    } else {
      fprintf(stderr,"alloc0_a: Error unable to allocate %ld bytes of space "
	      "for solvent_string\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  return(success);
}


