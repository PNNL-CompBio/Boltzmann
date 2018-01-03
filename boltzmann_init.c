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
#include "sort_compartments.h"
#include "unique_compartments.h"
#include "translate_compartments.h"
#include "sort_molecules.h"
#include "unique_molecules.h"
#include "print_molecules_dictionary.h"
#include "alloc3.h"
#include "set_compartment_ptrs.h"
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
  struct vgrng_state_struct *vgrng2_state;
  struct rxn_struct *reactions;
  struct istring_elem_struct *cur_molecules;
  struct istring_elem_struct *cur_cmpts;
  struct istring_elem_struct *cur_cmpt;
  char *cmpt_string;
  double *dg0s;
  double *free_energy;
  double *activities;
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
  int i;

  int nu_molecules;
  int padi;

  int oi;
  int ci;

  FILE *bndry_flux_fp;
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
      /*
	Open the log file.
      */
      state->lfp = fopen(state->log_file,"w");
      lfp = state->lfp;
      if (state->lfp == NULL) {
	fprintf(stderr,
		"boltzman_init unable to open log_file, %s, quitting.\n",
		state->log_file);
	fflush(stderr);
	success = 0;
      }
    }
  }
  if (success) {
    if (state->concs_out_file) {
      /*
	Open the concentrations output file.
      */
      state->concs_out_fp = fopen(state->concs_out_file,"w");
      if (state->concs_out_fp == NULL) {
	fprintf(stderr,
		"boltzman_init unable to open concs_out_file, %s, quitting.\n",
		state->concs_out_file);
	fflush(stderr);
	success = 0;
      }
    }
  }
  if (success) {
    if (state->rxn_lklhd_file) {
      /*
	Open the likelihoods output file.
      */
      state->rxn_lklhd_fp = fopen(state->rxn_lklhd_file,"w");
      if (state->rxn_lklhd_fp == NULL) {
	fprintf(stderr,
		"boltzman_init unable to open rxn_lklhd_file, %s, quitting.\n",
		state->rxn_lklhd_file);
	fflush(stderr);
	success = 0;
      }
    }
  }
  if (success) {
    if (state->free_energy_format > 0) {
      /*
	Open the free energy output file.
      */
      if (state->free_energy_file) {
	state->free_energy_fp = fopen(state->free_energy_file,"w");
	if (state->free_energy_fp == NULL) {
	  fprintf(stderr,
		  "boltzman_init unable to open free_energy_file, %s, quitting.\n",
		  state->free_energy_file);
	  fflush(stderr);
	  success = 0;
	}
      }
    }
  }
  if (success) {
    if (state->bndry_flux_file) {
      /*
	Open the boundary flux output file.
      */
      bndry_flux_fp = fopen(state->bndry_flux_file,"w");
      if (bndry_flux_fp == NULL) {
	fprintf(stderr,
		"boltzman_init unable to open bndry_flux_file, %s, quitting.\n",
		state->bndry_flux_file);
	fflush(stderr);
	success = 0;
      } else {
	state->bndry_flux_fp = bndry_flux_fp;
      }
    }
  }
  if (success) {
    /*
      Initialize the random number generator.
    */
    vgrng_state = state->vgrng_state;
    vgrng_start_steps = 1001;
    vgrng_start= vgrng_init(vgrng_state,vgrng_start_steps);
    vgrng2_state = state->vgrng2_state;
    vgrng_start_steps = 1042;
    vgrng_start= vgrng_init(vgrng2_state,vgrng_start_steps);
  }
  /*
    Allocate space for the reactions line buffer, and the rxn_file keywords.
  */
  if (success) {
    success = alloc1(state);
  }
  if (success) {
    /*
      Echo the input parameters to the log file.
    */
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
    First we need to sort the compartments.
  */
  if (success) {
    success = sort_compartments(&state->unsorted_cmpts,
				&state->sorted_cmpts,
				state->number_compartments);
  }
  if (success) {
    success = unique_compartments(state);
  }
  /*
    Now we need to assign the proper compartment numbers to the 
    unsorted molecules, using the compartment_indices field
    of the rxns_matrix.
  */
  if (success) {
    success = translate_compartments(state);
  }
  /*
    Now we need to sort the molecules.
  */
  if (success) {
    success = sort_molecules(&state->unsorted_molecules,
			    &state->sorted_molecules,
			    state->number_molecules);
  }
  if (success) {
    success = unique_molecules(state);
  }
  if (success) {
    nu_molecules = state->unique_molecules;
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
  /*
    Now in order to enable molecule/compartment lookup for
    read_initial_concentrations we need to set the compartment 
    pointers.
  */
  if (success) {
    success = set_compartment_ptrs(state);
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
  if (success) {
    fprintf(state->rxn_lklhd_fp,"iter\tentropy\tdg_forward\tforward_rxn_likelihood\treverse_rxn_likelihood\n");
    fprintf(state->rxn_lklhd_fp,"iter\tentropy\tdg_forward");
    reactions                   = state->reactions;
    for (i=0;i<state->number_reactions;i++) {
      fprintf(state->rxn_lklhd_fp,"\tf_%s\tr_%s",reactions->title,reactions->title);
      reactions += 1; /* Caution address arithmetic */
    }
    fprintf(state->rxn_lklhd_fp,"\n");
    fflush(state->rxn_lklhd_fp);
  }
  if (success) {
    if (state->free_energy_format > 0) {
      if (state->free_energy_format == 1) {
	fprintf(state->free_energy_fp,"negative_log_likelihoods\n");
      } else if (state->free_energy_format == 2) {
	fprintf(state->free_energy_fp,"free energy (KJ/mol)\n");
      } else if (state->free_energy_format == 3) {
	fprintf(state->free_energy_fp,"free energy (Kcal/mol)\n");
      }
      fprintf(state->free_energy_fp,"iter");
      reactions                   = state->reactions;
      for (i=0;i<state->number_reactions;i++) {
	fprintf(state->free_energy_fp,"\t%s",reactions->title);
	reactions += 1; /* Caution address arithmetic */
      }
      fprintf(state->free_energy_fp,"\n");
      fflush(state->free_energy_fp);
    }
  }
  /*
  if (success) {
    if (state->num_fixed_concs >0) {
      if (bndry_flux_fp) {
	cur_molecules = state->sorted_molecules;
	cur_cmpts     = state->sorted_cmpts;
	fprintf(bndry_flux_fp,"       iter       ");
	cmpt_string = NULL;
	oi          = -1;
	for (i=0;i<nu_molecules;i++) {
	  ci = cur_molecules->c_index;
	  if (ci != oi) {
	    oi = ci;
	    cur_cmpt = (struct istring_elem_struct *)&(cur_cmpts[ci]);
	    cmpt_string = cur_cmpt->string;
	  }
	  if (cur_molecules->variable == 0) {
	    if (ci != -1) {
	      fprintf(bndry_flux_fp,"\t%s:%s",
		      cur_molecules->string,cmpt_string);
	    } else {
	      fprintf(bndry_flux_fp,"\t%s",cur_molecules->string);
	    }
	  }
	  cur_molecules += 1; // Caution Address arithmetic.
	}
	fprintf(bndry_flux_fp,"\n");
      } // end if (bndry_flux_fp).
    
    
    } // end if (state->num_fixed_concs...) 
  }
  */
  if (success) {
    dg0s = state->dg0s;
    free_energy  = state->free_energy;
    for (i=0;i<state->number_reactions;i++) {
      free_energy[i] = dg0s[i];
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

  return(success);
}
