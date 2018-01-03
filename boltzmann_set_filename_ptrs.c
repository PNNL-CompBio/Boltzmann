/* boltzmann_set_filename_ptrs.c
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

void boltzmann_set_filename_ptrs(struct state_struct *statep) {
  /*
    Set the char * pointer fields of the state vector that
    point to filenames for input and output files and directories.
    Called by: alloc0,flatten_state.
  */

  int64_t max_file_name_len;
  max_file_name_len          = statep->max_filename_len;
  /*
    Caution Address arthmetic follows:
  */
  statep->reaction_file      = statep->params_file + max_file_name_len;
  statep->init_conc_file     = statep->reaction_file + max_file_name_len;
  statep->input_dir          = statep->init_conc_file + max_file_name_len;
  statep->output_file        = statep->input_dir + max_file_name_len;
  statep->log_file           = statep->output_file + max_file_name_len;
  statep->output_dir         = statep->log_file + max_file_name_len;
  statep->counts_out_file    = statep->output_dir + max_file_name_len;
  statep->rxn_lklhd_file     = statep->counts_out_file + max_file_name_len;
  statep->free_energy_file   = statep->rxn_lklhd_file + max_file_name_len;
  statep->restart_file       = statep->free_energy_file + max_file_name_len;
  statep->rxn_view_file      = statep->restart_file + max_file_name_len;
  statep->bndry_flux_file    = statep->rxn_view_file + max_file_name_len;
  statep->pseudoisomer_file  = statep->bndry_flux_file + max_file_name_len;
  statep->compartment_file   = statep->pseudoisomer_file + max_file_name_len;
  statep->sbml_file          = statep->compartment_file + max_file_name_len;
  statep->ms2js_file         = statep->sbml_file + max_file_name_len;
  statep->kg2js_file         = statep->ms2js_file + max_file_name_len;
  statep->rxn_echo_file      = statep->kg2js_file + max_file_name_len;
  statep->rxn_mat_file       = statep->rxn_echo_file + max_file_name_len;
  statep->dg0ke_file         = statep->rxn_mat_file + max_file_name_len;
  statep->dictionary_file    = statep->dg0ke_file + max_file_name_len;
}
