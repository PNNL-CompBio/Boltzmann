/* boltzmann_init.c
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

#include "djb_timing_b.h"
#include "boltzmann_structs.h"
#include "alloc0.h"
#include "read_params.h"
#include "vgrng_init.h"
#include "alloc1.h"
#include "echo_params.h"
#include "size_rxns_file.h"
#include "alloc2.h"
#include "parse_reactions_file.h"
#include "echo_reactions_file.h"
#include "sort_istrings.h"
#include "unique_molecules.h"
#include "print_molecules_dictionary.h"
#include "alloc3.h"
#include "read_initial_concentrations.h"
#include "form_molecules_matrix.h"
#include "compute_ke.h"
/*
#define DBG_BOLTZMANN_INIT  
*/


#include "boltzmann_init.h"
int boltzmann_init(char *param_file_name, struct state_struct **statep) {
  /*
    Initialize the reactions and data structures for boltzmann.
    Called by: boltzmann
    Calls:     alloc0,
               read_params,
	       vgrng_init,
	       alloc1,
	       echo_params,
	       size_rxns_file,
	       alloc2,
	       parse_reactions_file,
	       echo_reactions_file,
	       sort_istrings,
	       unique_molecules,
	       print_molecues_dictionary,
	       alloc3,
	       read_initial_concentrations,
	       form_molecules_matrix,
	       compute_ke
  */
  struct state_struct bltzs;
  struct state_struct *state;
  struct vgrng_state_struct *vgrng_state;
  int64_t align_len;
  int64_t align_mask;
  int64_t rxn_title_len;
  int64_t pathway_len;
  int64_t compartment_len;
  int64_t molecules_len;
  int64_t vgrng_start;
  int success;
  int num_state_files;

  int num_rxns;
  int num_molecules;

  int vgrng_start_steps;
  int padi;

  FILE *lfp;
  /*
    allocate space for the state struct.
  */
  success = alloc0(statep);
  if (success) {
    /*
      Read the input parameters file.
    */
    state = *statep;
    success = read_params(param_file_name,state);
  }
  if (success) {
    if (state->log_file) {
      state->lfp = fopen(state->log_file,"w");
      lfp = state->lfp;
      if (state->lfp == NULL) {
	fprintf(stderr,"boltzman unable to open log_file, %s, quitting.\n",
		state->log_file);
	fflush(stderr);
	success = 0;
      }
    }
  }
  if (success) {
    vgrng_state = state->vgrng_state;
    vgrng_start_steps = 1001;
    vgrng_start= vgrng_init(vgrng_state,vgrng_start_steps);
  }
  /*
    Allocate space for the reactions line buffer, and the rxn_file keywords.
  */
  if (success) {
    success = alloc1(state);
  }
  if (success) {
    success = echo_params(state->lfp,state);
  }
  /*
    Read reactions file to count molecules and reactions
    Then allocate space then read for real.
  */
  if (success) {
    success = size_rxns_file(state);
  }
  if (success) {
    lfp = state->lfp;
#ifdef DBG_BOLTZMANN_INIT
    if (lfp) {
      fprintf(lfp,"boltzmann_init: after startup, success = %d\n",success);
      fprintf(lfp,"boltzmann_init: rxns = %d, cmpts = %d, molecules = %d, "
	      "reaction_file_length = %ld\n",
	      state->number_reactions,state->number_compartments,
	      state->number_molecules,
	      state->reaction_file_length);
      fprintf(lfp,"boltzmann_init: molecules_len = %ld, reaction_title_len = %ld, "
	      "pathway_len = %ld, compartment_len = %ld\n",
	      state->molecules_len,
	      state->rxn_title_len,
	      state->pathway_len,
	      state->compartment_len);
      fflush(lfp);
    }
#endif
  }
  /*
    At this point in time we can compute how much space the aligned
    reaction titles, pathway descriptions, compartments and molecules
    verbage will take.
    We want an uppercase version of the molecules, so we need two 
    molecules copies. Also we have an upperbound of the number molecules,
    state->number_molecules and we can allocate for the molecules sorting.
  */
  if (success) {
    success = alloc2(state);
  }
  if (success) {
    success = parse_reactions_file(state);
  }
  if (success) {
    success = echo_reactions_file(state);
  }
  /*
    Now we need to sort the molecules.
  */
  if (success) {
    success = sort_istrings(&state->unsorted_molecules,
			    &state->sorted_molecules,
			    state->number_molecules);
  }
  if (success) {
    success = unique_molecules(state);
  }
  if (success) {
    success = print_molecules_dictionary(state);
  }
  /*
    Now we need to allocate space for the concentrations,
    transpose the reactions matrix to the molecules matrix,
    and read in the intial concentrations.
  */
  if (success) {
    success = alloc3(state);
  }
  if (success) {
    success = read_initial_concentrations(state);
  }
  if (success) {
    success = form_molecules_matrix(state);
  }
  if (success) {
    success = compute_ke(state);
  }
  return(success);
}
