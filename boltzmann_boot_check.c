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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <float.h>
#include <signal.h>
#include <unistd.h>
#include <errno.h>
#include <sys/mman.h>


#include "djb_timing_b.h"
#include "boltzmann_structs.h"
#include "boltzmann_rep_state_i.h"
#include "echo_reactions_file.h"
/*
#define DBG_BOLTZMANN_BOOT_CHECK 1  
*/
#include "boltzmann_boot_check.h"
int boltzmann_boot_check(int64_t *super_statep, FILE *lfp) {
  struct  state_struct *lstate;
  double  *concs;
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
  int64_t molecule_map_starts_offset;
  int64_t molecule_map_offset;
  int64_t compartment_map_offset;
  int64_t molecule_names_offset;
  int64_t compartment_names_offset;
  int64_t molecules_text_offset;
  int64_t compartments_text_offset;
  int64_t map_index;
  int64_t molecule;
  int64_t compartment;
  int64_t one_l;
  
  int     success;
  int     j;
  success = 1;
  one_l   = (int64_t)1;
  super_statec = (char*)super_statep;
  super_statel = (int64_t*)super_statep;
  num_reaction_files = super_statel[0];
  state_offsets_sizes = &super_statel[super_statel[12]];
  molecule_map_starts_offset  = super_statel[13];
  molecule_map_offset      = super_statel[14];
  molecule_names_offset    = super_statel[15];
  compartment_map_offset   = super_statel[16];
  compartment_names_offset = super_statel[17];
  molecules_text_offset    = super_statel[18];
  compartments_text_offset = super_statel[19];
  molecule_map_starts      = &super_statel[molecule_map_starts_offset];
  molecule_map             = &super_statel[molecule_map_offset];
  compartment_map          = &super_statel[compartment_map_offset];
  molecule_names           = &super_statel[molecule_names_offset];
  compartment_names        = &super_statel[compartment_names_offset];
  molecules_text           = &super_statec[molecules_text_offset];
  compartments_text        = &super_statec[compartments_text_offset];
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
      success = boltzmann_rep_state_i(super_statep,(int)i,&lstate);
      if (success) {
	/*
	  Print out the reactions in the i'th reaction file from
	  its state structre.
	*/
	fprintf (lfp,"reaction file: %s, state size = %ld\n",
		 lstate->reaction_file,lsize);
	success = echo_reactions_file(lstate,lfp);
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
	concs = lstate->current_concentrations;
	fprintf (lfp,"\n\n Initial concentrations from %s\n",
		 lstate->init_conc_file);
	for (j=0;j<lstate->nunique_molecules;j++) {
	  molecule    = *local_molecules;
	  compartment = *local_compartments;
	  if (compartment > 0) {
	    fprintf(lfp,"%s:%s\t%le\n",
		    (char*)&molecules_text[molecule_names[molecule]],
		    (char*)&compartments_text[compartment_names[compartment]],
		    concs[j]);
	  } else {
	    fprintf(lfp,"%s\t%le\n",
		    (char*)&molecules_text[molecule_names[molecule]],
		    concs[j]);
	  }
	  local_molecules += 1; /* Caution address arithmetic */
	  local_compartments += 1; /* Caution address arithmetic */
	}
	fprintf(lfp,"\n---------------------------------\n");
      }
      free(lstate);
    }
  }
  return(success);
}
      
