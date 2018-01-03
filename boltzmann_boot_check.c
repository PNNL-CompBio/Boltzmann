/* boltzmann_boot_check.c
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
#include "flatten_super_state.h"
#include "boltzmann_rep_state_i.h"
#include "echo_reactions_file.h"
/*
#define DBG_BOLTZMANN_BOOT_CHECK 1  
*/
#include "boltzmann_boot_check.h"
int boltzmann_boot_check(struct super_state_struct *super_statep, FILE *lfp) {
  struct  super_state_pointers_struct ssps;
  struct  super_state_pointers_struct *super_state_pointers;
  struct  state_struct *lstate;
  double  *counts;
  int64_t *super_statel;
  int64_t *state_offsets_sizes;
  int64_t *molecule_map_starts;
  int64_t *molecule_map;
  int64_t *compartment_map;
  int64_t *molecule_names;
  int64_t *compartment_names;
  int64_t *local_molecules;
  int64_t *local_compartments;
  char    *super_statec;
  char    *molecules_text;
  char    *compartments_text;
  int64_t num_reaction_files;
  int64_t i;
  int64_t loffset;
  int64_t lsize;
  int64_t map_index;
  int64_t molecule;
  int64_t compartment;
  int64_t one_l;
  int64_t ask_for;
  
  int     success;
  int     j;
  success = 1;
  one_l   = (int64_t)1;
  super_statec = (char*)super_statep;
  super_statel = (int64_t*)super_statep;
  num_reaction_files = super_statep->number_of_reaction_files;
  ask_for = sizeof(ssps);
  super_state_pointers = (struct super_state_pointers_struct *)calloc(one_l,ask_for);
  if (super_state_pointers != NULL) {
    success = flatten_super_state(super_statep,super_state_pointers);
  } else {
    fprintf(stderr,"boltzmann_boot_check: Error, unable to allocate %ld "
	    "bytes for super_state_pointers\n",ask_for);
    fflush(stderr);
    success = 0;
  }
  if (success) {
    state_offsets_sizes      = super_state_pointers->state_offsets_sizes;
    molecule_map_starts      = super_state_pointers->molecule_map_starts;
    molecule_map             = super_state_pointers->molecule_map;
    compartment_map          = super_state_pointers->compartment_map;
    molecule_names           = super_state_pointers->molecule_names;
    compartment_names        = super_state_pointers->compartment_names;
    molecules_text           = super_state_pointers->molecules_text;
    compartments_text        = super_state_pointers->compartments_text;
    local_molecules          = &molecule_map[0];
    local_compartments       = &compartment_map[0];
    if (lfp) {
      for (i=0;i<num_reaction_files;i++) {
	/*
	  Allocate space for a copy of the state_struct for the i'th
	  reaction file, returning pointer to it in lstate, and
	  initialize it and set its self pointers for use in
	  printing out concentrations.
	*/
	lstate = NULL;
	success = boltzmann_rep_state_i(super_statep,(int)i,&lstate);
	if (success) {
	  /*
	    Print out the reactions in the i'th reaction file from
	    its state structre.
	  */
	  fprintf (lfp,"reaction file: %s, state size = %ld\n",
		   lstate->reaction_file,lsize);
	  lstate->lfp = lfp;
	  success = echo_reactions_file(lstate);
	}
	/*
	  Now print out the intial concentrations using the
	  global molecules_map  and compartments_map vectors
	  to translate local molecules and compartment numbers
	  to the global dictionary.
	*/
	if (success) {
	  /*
	    map_index = map_indices[i];
	    local_molecules = &molecule_map[map_index];
	    local_compartments = &compartment_map[map_index];
	  */
	  counts = lstate->current_counts;
	  fprintf (lfp,"\n\n Initial counts from %s\n",
		   lstate->init_conc_file);
	  for (j=0;j<lstate->nunique_molecules;j++) {
	    molecule    = *local_molecules;
	    compartment = *local_compartments;
	    if (compartment > 0) {
	      fprintf(lfp,"%s:%s\t%le\n",
		      (char*)&molecules_text[molecule_names[molecule]],
		      (char*)&compartments_text[compartment_names[compartment]],
		      counts[j]);
	    } else {
	      fprintf(lfp,"%s\t%le\n",
		      (char*)&molecules_text[molecule_names[molecule]],
		      counts[j]);
	    }
	    local_molecules += 1; /* Caution address arithmetic */
	    local_compartments += 1; /* Caution address arithmetic */
	  }
	  fprintf(lfp,"\n---------------------------------\n");
	}
	free(lstate);
      } /* end for(i...) */
    }
  }
  return(success);
}
      
