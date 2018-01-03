/* rxn_map_init.c
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
#include "alloc0.h"
#include "read_params.h"
#include "open_output_files.h"
#include "echo_params.h"
#include "size_rxns_file.h"
#include "alloc2.h"
#include "parse_reactions_file.h"
#include "echo_reactions_file.h"
#include "sort_compartments.h"
#include "unique_compartments.h"
#include "translate_compartments.h"
#include "sort_molecules.h"
#include "unique_molecules.h"
#include "print_molecules_dictionary.h"
#include "alloc3.h"
#include "alloc4.h"
#include "set_compartment_ptrs.h"
#include "read_initial_concentrations.h"
#include "form_molecules_matrix.h"
/*
#define DBG_RXN_MAP_INIT 1
*/
#include "boltzmann_init.h"
int rxn_map_init(char *param_file_name, struct state_struct **statep,
		 struct molecules_matrix_struct **molecules_matrix_p) {
  /*
    Initialize the reactions and data structures for rxn_map.
    Called by: rxn_map
    Calls:     alloc0,
               read_params,
	       open_output_files,
	       echo_params,
	       size_rxns_file,
	       alloc2,
	       parse_reactions_file,
	       echo_reactions_file,
	       sort_compartments,
	       unique_compartments,
	       translate_compartments,
	       sort_molecules,
	       unique_molecules,
	       print_molecules_dictionary,
	       alloc3,
	       read_initial_concentrations,
	       form_molecules_matrix,
	       print_rxn_likelihoods_header
	       print_free_energy_header
  */
  struct state_struct *state;
  struct molecules_matrix_struct *molecules_matrix;
  int64_t *transpose_work;
  int success;
  int print_output;

  FILE *start_stop_fp;
  FILE *rxn_echo_fp;
#ifdef DBG_BOLTZMANN_INIT
  FILE *lfp;
  FILE *efp;
#endif
  /*
    allocate space for the state struct.
    Allocate space for the reactions line buffer, and the rxn_file keywords.
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
    print_output = state->print_output;
    if (state->log_file) {
      /*
	Open the log file.
      */
      state->lfp = fopen(state->log_file,"w");
      if (state->lfp == NULL) {
	fprintf(stderr,
		"rxn_map_init: unable to open log_file, %s, quitting.\n",
		state->log_file);
	fflush(stderr);
	success = 0;
      }
    }
  }
  if (success) {
    /*
      Echo the input parameters to the log file.
    */
    if (print_output) {
      success = echo_params(state->lfp,state);
    }
  }
  /*
    Read reactions file to count molecules and reactions
    Then allocate space then read for real.
  */
  if (success) {
    success = size_rxns_file(state,
			     state->reaction_file);
  }
  if (success) {
#ifdef DBG_BOLTZMANN_INIT
    lfp = state->lfp;
    if (lfp) {
      fprintf(lfp,"boltzmann_init: after startup, success = %d\n",success);
      fprintf(lfp,"boltzmann_init: rxns = %ld, cmpts = %ld, molecules = %ld, "
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
    molecules copies. Also we have an upperbound of the number of molecules,
    state->number_molecules, and we can allocate space for the molecules 
    sorting.
  */
  if (success) {
    success = alloc2(state);
  }
  if (success) {
    /*
      Read and parse the reactions file.
    */
    success = parse_reactions_file(state,state->reaction_file);
  }
  if (success) {
    /*
      Echo the reactions to the log file.
    */
    if (print_output) {
      rxn_echo_fp = fopen("rxns.echo","w+");
      if (rxn_echo_fp == NULL) {
	fprintf(stderr,
		"rxn_map_init: Error could not open rxns.echo file.\n");
	success = 0;
      }
      if (success) {
	success = echo_reactions_file(state);
	fclose(rxn_echo_fp);
      }
    }
  }
  /*
    First we need to sort the compartments.
  */
  if (success) {
    success = sort_compartments(state->unsorted_cmpts,
				state->sorted_compartments,
				state->compartment_text,
				state->number_compartments);
  }
  /*
    Then we extract the unique compartments.
    and fill the compartment_indices vector of the rxns_matrix structure.
  */
  if (success) {
    success = unique_compartments(state);
  }
  /*
    Now we need to assign the proper compartment numbers to the 
    unsorted molecules, using the compartment_indices field
    of the rxns_matrix structure.
  */
  if (success) {
    success = translate_compartments(state);
  }
  /*
    Now we need to sort the molecules, by compartment and name.
  */
  if (success) {
    success = sort_molecules(state->unsorted_molecules,
			     state->sorted_molecules,
			     state->molecules_text,
			     state->number_molecules);
  }
  /*
    Then we extract the unique molecules and set the
    molecules_indices field of the rxn_matrix struct.
  */
  if (success) {
    success = unique_molecules(state);
  }
  /*
    Print the molecules dictionary and the header lines for 
    the concentrations output file.
  */
  if (success) {
    if (print_output) {
      success = print_molecules_dictionary(state);
    }
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
    success = alloc4(state,&molecules_matrix,&transpose_work);
  }
  /*
    Now in order to enable molecule/compartment lookup for
    read_initial_concentrations we need to set the compartment 
    pointers, these are pointers into the list of sorted molecules
    All the molecules within a compartment are adjacent in the
    sorted molecules list. This call just sets pointers to the
    first molecule in each compartment, and one past the last
    molecule so that molecules in compartment i in the sorted molecules list
    are in positions compartment_ptrs[i]:compartment_ptrs[i+1]-1 
    inclusive.
  */
  if (success) {
    success = set_compartment_ptrs(state);
  }
  /*
    Open the  "initial concentrations file" which will really be
    a set of lines with start_molecule:compartment end_molecule:compartment
    lines.
  */
  if (success) {
    start_stop_fp = fopen(state->init_conc_file,"r");
    if (start_stop_fp == NULL) {
      fprintf(stderr,"rxn_map_init: Error unable to open file "
	     "state->init_conc_file,%s\n",state->init_conc_file);
      fflush(stderr);
      success = 0;
    } else {
      state->conc_fp = start_stop_fp;
    }
  }
  /*
    Compute the molecules matrix.
  */
  if (success) {
    success = form_molecules_matrix(state,molecules_matrix,transpose_work);
    *molecules_matrix_p = molecules_matrix;
  }
  return(success);
}
