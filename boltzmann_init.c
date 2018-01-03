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
#include "boltzmann_structs.h"
#include "alloc0.h"
#include "read_params.h"
#include "open_output_files.h"
#include "vgrng_init.h"
#include "echo_params.h"
#include "size_rxns_file.h"
#include "size_pseudoisomer_file.h"
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
#include "set_compartment_ptrs.h"
#include "read_initial_concentrations.h"
#include "compute_standard_energies.h"
#include "compute_ke.h"
#include "compute_kss.h"
#include "print_dg0_ke.h"
#include "zero_solvent_coefficients.h"
#include "print_rxn_likelihoods_header.h"
#include "print_free_energy_header.h"
#include "flatten_state.h"
#include "print_reactions_matrix.h"
#include "free_boot_state.h"
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
	       open_output_files,
	       vgrng_init,
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
	       compute_standard_energies,
	       compute_ke,
	       compute_kss,
	       print_rxn_likelihoods_header,
	       print_free_energy_header
  */
  struct state_struct *boot_state;
  struct state_struct *state;
  struct state_struct *stateq;
  struct formation_energy_struct *formation_energies;
  struct vgrng_state_struct *vgrng_state;
  struct vgrng_state_struct *vgrng2_state;
  double *dg0s;
  double *free_energy;
  double *activities;
  int64_t vgrng_start;
  int64_t i;

  int success;
  int vgrng_start_steps;

  int print_output;
  int padi;
  
  FILE *rxn_echo_fp;
  FILE *bndry_flux_fp;
  FILE *lfp;
  /*
    allocate space for the state struct.
    Allocate space for the reactions line buffer, and the rxn_file keywords.
  */
  success = alloc0(&boot_state);
  if (success) {
    /*
      Read the input parameters file.
    */
    state = boot_state;
    success = read_params(param_file_name,state);
  }
  if (success) {
    print_output = (int)state->print_output;
    if (print_output) {
      success = open_output_files(state);
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
  if (success) {
    success = size_rxns_file(state,state->reaction_file);
  }
  if (success) {
    success = alloc2(state);
  }
  /*
    Read reactions file to count molecules and reactions
    Then allocate space then read for real.
  */
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
    /*
      Read and parse the reactions file.
    */
    success = parse_reactions_file(state,state->reaction_file);
  }
  /*
    First we need to sort the compartments.
  */
  if (success) {
    success = sort_compartments(state->unsorted_cmpts,
				state->sorted_cmpts,
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
    Also set the solvent_pos field of state.
  */
  if (success) {
    success = unique_molecules(state);
  }
  /*
    Now we need to allocate space for the counts, concentrations,
    and read in the intial concentrations converting them to counts.
  */
  if (success) {
    success = alloc3(state);
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
    Read initial concentrations.
    And print them to the counts output file.
  */
  if (success) {
    success = read_initial_concentrations(state);
  }
  /*
    Compute the molecules matrix.
  if (success) {
    success = form_molecules_matrix(state);
  }
  */
  /*
    Compute the reaction energies of formation if called for.
  */
  if (success) {
    if (state->use_pseudoisomers) {
      formation_energies = NULL;
      success = compute_standard_energies(state,&formation_energies);
    }
  }
  if (success) {
    if (print_output) {
      /*
	Echo the reactions to the log file.
      */
      rxn_echo_fp = fopen("rxns.echo","w+");
      if (rxn_echo_fp == NULL) {
	fprintf(stderr,
		"echo_reactions_file: Error could not open rxns.echo file.\n");
	success = 0;
      }
      if (success) {
	success = echo_reactions_file(state,rxn_echo_fp);
	fclose(rxn_echo_fp);
      }

    }
  }
  /*
    Print the molecules dictionary and the header lines for 
    the counts output file.
  */
  if (success) {
    if (print_output) {
      success = print_molecules_dictionary(state);
    }
  }
  /*
    Compute the reaction ke's.
  */
  if (success) {
    success = compute_ke(state);
  }
  /*
    Compute the reaction kss's.
  */
  if (success) {
    success = compute_kss(state);
  }
  if (success) {
    if (print_output) {
      success = print_dg0_ke(state);
    }
  }
  /*
    At this juncture we have echoed the reactions file if requested and
    need to zero out the coefficients in the reaction matrix that
    correspond to the solvent molecule (by default H2O) so as not to
    have it influence the computation of likelihoods, nor change
    concentration (see rxn_likelihood.c and comment in rxn_conc_update.c)
  */
  if (success) {
    success = zero_solvent_coefficients(state);
  }
  if (success) {
    /*
      Print the header lines for the reaction likelihoods output file.
    */
    if (print_output) {
      print_rxn_likelihoods_header(state);
    }
  }
  if (success) {
    if (state->free_energy_format > (int64_t)0) {
      if (print_output) {
	/*
	  Print the header lines for the free energy output file.
	*/
	print_free_energy_header(state);
      }
    }
  }
  /*
    If use_activities has not been turned on, set all activities to
    1.0 so that all reactions are fully active.
  */
  if (success) {
    activities = state->activities;
    if (state->use_activities == 0) {
      for (i=0;i<state->number_reactions;i++) {
	activities[i] = 1.0;
      }
    }
  }
  if (success) {
    /*
      Initialize the random number generators,
      setting the vgrng_state and vgrng2_state fields of
      the state structure.
    */
    vgrng_state = state->vgrng_state;
    vgrng_start_steps = 1001;
    vgrng_start= vgrng_init(vgrng_state,vgrng_start_steps);
    vgrng2_state = state->vgrng2_state;
    vgrng_start_steps = 1042;
    vgrng_start= vgrng_init(vgrng2_state,vgrng_start_steps);
  }
  if (success) {
    if (boot_state->print_output) {
      boot_state->rxn_view_hist_length = ((int64_t)(boot_state->record_steps + boot_state->rxn_view_freq -1)/boot_state->rxn_view_freq) + (int64_t)1;
    } else {
      boot_state->rxn_view_hist_length = 0;
    }
    /*
    *statep = boot_state;
    */
    stateq  = NULL;
    boot_state->workspace_base = NULL;
    success = flatten_state(boot_state,&stateq);
    if (success) {
      if (state->print_output) {
	success = print_reactions_matrix(stateq);
      }
    }
    *statep = stateq;
    if (success) {
      success = free_boot_state(&boot_state);
    }
  }
  return(success);
}
