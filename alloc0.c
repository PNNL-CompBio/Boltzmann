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

#include "alloc0.h"
int alloc0(struct state_struct **statep) {
  /*
    Allocate space for boltzmann_state data and the associated state
    file names. And the parameter file input buffer line, and
    the random number generator state.
    Called by: boltzmann_init
    Calls    : calloc, fprintf, fflush (intrinsic)
  */
  struct state_struct bltzs;
  struct state_struct *state;
  int64_t *rxn_file_keyword_lengths;
  char    *rxn_keyword_buff;
  char    *solvent_string;
  char    **rxn_keywords;
  int64_t max_file_name_len;
  int64_t max_param_line_len;
  int64_t ask_for;
  int64_t one_l;
  int64_t usage;
  int64_t num_state_files;

  int success;

  int num_rxn_file_keywords;
  int keyword_buffer_length;

  ask_for           = (int64_t)sizeof(bltzs);
  one_l             = (int64_t)1;
  max_file_name_len = (int64_t)4096;
  max_param_line_len = (int64_t)4096;
  success           = 1;
  usage             = ask_for;
  state             = (struct state_struct *)calloc(one_l,ask_for);
  *statep           = state;
  if (state == NULL) {
    success = 0;
    fprintf(stderr,
	    "boltzmann: unable to allocate %lld bytes for state structure.\n",
	    ask_for);
    fflush(stderr);
  }
  if (success) {
    /*
      Right now we only have 28 file_names 
      but we leave space for additional ones.
    */
    state->num_files = (int64_t)32;
    num_state_files   = state->num_files;
    state->max_filename_len = max_file_name_len;
    ask_for = ((int64_t)num_state_files) * max_file_name_len;
    usage   += ask_for;
    state->params_file = (char *)calloc(one_l,ask_for);
    if (state->params_file == NULL) {
      success = 0;
      fprintf(stderr,
	      "alloc0: unable to allocate %lld bytes for state file names.\n",
	      ask_for);
      fflush(stderr);
    }
  }
  if (success) {
    state->max_filename_len   = max_file_name_len;
    boltzmann_set_filename_ptrs(state);
    /*
    state->reaction_file      = state->params_file + max_file_name_len;
    state->init_conc_file     = state->reaction_file + max_file_name_len;
    state->input_dir          = state->init_conc_file + max_file_name_len;
    state->output_file        = state->input_dir + max_file_name_len;
    state->log_file           = state->output_file + max_file_name_len;
    state->output_dir         = state->log_file + max_file_name_len;
    state->counts_out_file    = state->output_dir + max_file_name_len;
    state->ode_concs_file     = state->counts_out_file + max_file_name_len;
    state->rxn_lklhd_file     = state->ode_concs_file + max_file_name_len;
    state->free_energy_file   = state->rxn_lklhd_file + max_file_name_len;
    state->restart_file       = state->free_energy_file + max_file_name_len;
    state->rxn_view_file      = state->restart_file + max_file_name_len;
    state->bndry_flux_file    = state->rxn_view_file + max_file_name_len;
    state->pseudoisomer_file  = state->bndry_flux_file + max_file_name_len;
    state->compartment_file   = state->pseudoisomer_file + max_file_name_len;
    state->sbml_file          = state->compartment_file + max_file_name_len;
    state->ms2js_file         = state->sbml_file + max_file_name_len;
    state->kg2js_file         = state->ms2js_file + max_file_name_len;
    */

    state->max_param_line_len = max_param_line_len;
    ask_for                    = max_param_line_len << 1;
    usage                      += ask_for;
    state->param_buffer       = (char *)calloc(one_l,ask_for);
    if (state->param_buffer == NULL) {
      success = 0;
      fprintf(stderr,
	      "alloc0: unable to allocate %lld bytes for state->param_buffer.\n",
	      ask_for);
      fflush(stderr);
    } 
  }
  if (success) {
    ask_for = ((int64_t)64) * sizeof(char);
    usage += ask_for;
    solvent_string = (char *)calloc(one_l,ask_for);
    if (solvent_string) {
      state->solvent_string = solvent_string;
    } else {
      fprintf(stderr,"alloc0: Error unable to allocate %lld bytes of space "
	      "for solvent_string\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  /*
    Allocate space for processing the reactions file.
  */
  if (success) {
    keyword_buffer_length = (int64_t)256;
    state->keyword_buffer_length = keyword_buffer_length;
    ask_for = keyword_buffer_length;
    usage   += ask_for;
    rxn_keyword_buff = (char *)calloc(one_l,ask_for);
    if (rxn_keyword_buff) {
      state->rxn_file_keyword_buffer = rxn_keyword_buff;
    } else {
      fprintf(stderr,
	      "alloc0: Error, unable to allocate %lld bytes for "
	      "rxn_file_keyword_buffer\n",
	      ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  if (success) {
    num_rxn_file_keywords = 14;
    state->num_rxn_file_keywords = num_rxn_file_keywords;
    ask_for = ((int64_t)num_rxn_file_keywords) * ((int64_t)sizeof(char *));
    usage   += ask_for;
    rxn_keywords = (char **)calloc(one_l,ask_for);
    if (rxn_keywords) {
      state->rxn_file_keywords = rxn_keywords;
    } else {
      fprintf(stderr,
	      "alloc0: Error, unable to allocate %lld bytes for "
	      "rxn_file_keywords\n",
	      ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  if (success) {
    ask_for = ((int64_t)num_rxn_file_keywords) * ((int64_t)sizeof(int64_t));
    usage += ask_for;
    rxn_file_keyword_lengths = (int64_t *)calloc(one_l,ask_for);
    if (rxn_file_keyword_lengths) {
      state->rxn_file_keyword_lengths = rxn_file_keyword_lengths;
    } else {
      fprintf(stderr,
	      "alloc0: Error, unable to allocate %lld bytes for "
	      "rxn_file_keyword_lengths\n",
	      ask_for);
      fflush(stderr);
      success = 0;
    }
    state->usage = usage;
  }
  return(success);
}
