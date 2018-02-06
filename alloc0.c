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
#include "boltzmann_cvodes_headers.h"
#include "cvodes_params_struct.h"
#include "alloc0_a.h"

#include "alloc0.h"
int alloc0(struct state_struct **statep, int setup) {
  /*
    Allocate space for boltzmann_state data and the associated state
    file names. And the parameter file input buffer line, and
    the random number generator state.
    If setup is 0 only allocate the state struct, as the
    other pieces are only used in setp. 
    If setup is 1 allocate all the pieces
    Called by: boltzmann_init, boltzmann_flatten_state
    Calls    : calloc, fprintf, fflush (intrinsic)
    Sets the following fields in state:

      set in alloc0_a:
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
	ode_counts_file,

        solvent_string

      Set after alloc0_a call

      version_no,
      max_param_line_len,
      param_buffer,
      vgrng_state,
      vgnrg2_state,
      num_rxn_file_keywords,
      rxn_file_keywords,
      keyword_buffer_length,
      rxn_file_keyword_buffer,
      rxn_file_keyword_lengths
  */
  struct state_struct bltzs;
  struct state_struct *state;
  struct vgrng_state_struct vss;
  struct cvodes_params_struct cps;
  struct ode23tb_params_struct ops;
  int64_t *rxn_file_keyword_lengths;
  char    *rxn_keyword_buff;
  char    **rxn_keywords;
  int64_t max_file_name_len;
  int64_t max_param_line_len;
  int64_t ask_for;
  int64_t one_l;
  int64_t usage;

  int success;

  int num_rxn_file_keywords;
  int keyword_buffer_length;

  ask_for           = (int64_t)sizeof(bltzs);
  one_l             = (int64_t)1;
  max_file_name_len = (int64_t)128;
  max_param_line_len = (int64_t)128;
  success           = 1;
  usage             = ask_for;
  state             = (struct state_struct *)calloc(one_l,ask_for);
  *statep           = state;
  if (state == NULL) {
    success = 0;
    fprintf(stderr,
	    "boltzmann: unable to allocate %ld bytes for state structure.\n",
	    ask_for);
    fflush(stderr);
  }
  if (success) {
    state->usage = usage;
    /*
      Allocate space for filename strings and solvent string
      These are considered auxilliary strings, needed only for setup
      and printing, not for running.
      This sets state->max_filename_len, the num_files parameter
      and allocates space for the filenames and the solvent string.
    */
    success = alloc0_a(state);
    usage = state->usage;
  }
  if (success && setup) {
    state->version_no = 5080; /* Checked in revision of state_struct.h */
    /*
      Allocate_space for reading the parameter file.
    */
    if (success) {
      state->max_param_line_len = max_param_line_len;
      ask_for                    = max_param_line_len << 1;
      usage                      += ask_for;
      state->param_buffer       = (char *)calloc(one_l,ask_for);
      if (state->param_buffer == NULL) {
	success = 0;
	fprintf(stderr,
		"alloc0: unable to allocate %ld bytes for state->param_buffer.\n",
		ask_for);
	fflush(stderr);
      } 
    }
    /*
      Allocate space for the random number generator states, vgrng_state
      and vgrng2_state.
    */
    if (success) {
      ask_for = (int64_t)sizeof(vss);
      usage += ask_for;
      state->vgrng_state = (struct vgrng_state_struct *)calloc(one_l,ask_for);
      if (state->vgrng_state == NULL) {
	success = 0;
	fprintf(stderr,
		"alloc0: unable to allocate %ld bytes for state->vgrng_state.\n",
		ask_for);
	fflush(stderr);
      }
    }
    if (success) {
      ask_for = (int64_t)sizeof(vss);
      usage += ask_for;
      state->vgrng2_state = (struct vgrng_state_struct *)calloc(one_l,ask_for);
      if (state->vgrng2_state == NULL) {
	success = 0;
	fprintf(stderr,
		"alloc0: unable to allocate %ld bytes for state->vgrng2_state.\n",
		ask_for);
	fflush(stderr);
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
		"alloc0: Error, unable to allocate %ld bytes for "
		"rxn_file_keyword_buffer\n",
		ask_for);
	fflush(stderr);
	success = 0;
      }
    }
    if (success) {
      num_rxn_file_keywords = 16;
      state->num_rxn_file_keywords = num_rxn_file_keywords;
      ask_for = ((int64_t)num_rxn_file_keywords) * ((int64_t)sizeof(char *));
      usage   += ask_for;
      rxn_keywords = (char **)calloc(one_l,ask_for);
      if (rxn_keywords) {
	state->rxn_file_keywords = rxn_keywords;
      } else {
	fprintf(stderr,
		"alloc0: Error, unable to allocate %ld bytes for "
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
	state->rxn_file_keyword_lengths = (int64_t *)rxn_file_keyword_lengths;
      } else {
	fprintf(stderr,
		"alloc0: Error, unable to allocate %ld bytes for "
		"rxn_file_keyword_lengths\n",
		ask_for);
	fflush(stderr);
	success = 0;
      }
    }
    /*
      Allocate space for the cvodes params struct.
    */
    if (success) {
      ask_for = (int64_t)sizeof(cps);
      usage+= ask_for;
      state->cvodes_params = (struct cvodes_params_struct *)calloc(one_l,ask_for);
      if (state->cvodes_params == NULL) {
	fprintf(stderr,"alloc0: Error, unalbe to allocate %ld bytes for "
		"cvodes_params\n",ask_for);
	fflush(stderr);
	success = 0;
      } else {
	state->cvodes_params_size = ask_for;
      }
    }
    /*
      Allocate space for the ode23tb params struct.
    */
    if (success) {
      ask_for = (int64_t)sizeof(ops);
      usage += ask_for;
      state->ode23tb_params = (struct ode23tb_params_sruct *)calloc(one_l,ask_for);
      if (state->ode23tb_params == NULL) {
	fprintf(stderr,"alloc0: Error, unalbe to allocate %ld bytes for "
		"ode23tb_params\n",ask_for);
	fflush(stderr);
	success = 0;
      } else {
	state->ode23tb_params_size = ask_for;
      }
    }
  } /* end if setup */
  state->usage = usage;
  return(success);
}
