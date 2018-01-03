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
  max_file_name_len = (int64_t)128;
  max_param_line_len = (int64_t)128;
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
  if (setup) {
    /*
      Allocate space for filename strings and solvent string
      These are considered auxilliary strings, needed only for setup
      and printing, not for running.
      This sets state->max_filename_len, the num_files parameter
      and allocates space for the filenames and the solvent string.
    */
    if (success) {
      state->usage = usage;
      success = alloc0_a(state);
      usage = state->usage;
    }
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
		"alloc0: unable to allocate %lld bytes for state->param_buffer.\n",
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
		"alloc0: Error, unable to allocate %lld bytes for "
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
    }
  } /* end if setup */
  state->usage = usage;
  return(success);
}
