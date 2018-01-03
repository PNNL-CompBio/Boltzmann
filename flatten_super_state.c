/* flatten_super_state.c
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
#include <unistd.h>
#include <string.h>
#include <strings.h>
#include <float.h>
#include <signal.h>

#include "boltzmann_structs.h"

#include "flatten_super_state.h"
/*
  Set the pointer fields of the super_state_struct.
  Called by: boltzmann_boot.
*/
int flatten_super_state (struct super_state_struct *super_statep,
			 struct super_state_pointers_struct *super_state_pointers) {
  int64_t *super_statel;
  char    *super_statec;
  int64_t state_offsets_sizes_offset;
  int64_t molecule_map_starts_offset;
  int64_t molecule_map_offset;
  int64_t molecule_names_offset;
  int64_t compartment_map_offset;
  int64_t compartment_names_offset;
  int64_t molecules_text_offset;
  int64_t compartments_text_offset;
  int success;
  int padi;
  success = 1;
  super_statec = (char*)super_statep;
  super_statel = (int64_t*)super_statep;
  state_offsets_sizes_offset  = super_statep->state_offsets_sizes_offset_in_words;  
  molecule_map_starts_offset  = super_statep->molecule_map_starts_offset_in_words;
  molecule_map_offset         = super_statep->molecule_map_offset_in_words;
  molecule_names_offset       = super_statep->molecule_names_offset_in_words;
  compartment_map_offset      = super_statep->compartment_map_offset_in_words;
  compartment_names_offset    = super_statep->compartment_names_offset_in_words;
  molecules_text_offset        = super_statep->molecules_text_offset_in_bytes;
  compartments_text_offset     = super_statep->compartments_text_offset_in_bytes;
  super_state_pointers->state_offsets_sizes = &super_statel[state_offsets_sizes_offset];
  super_state_pointers->molecule_map_starts = &super_statel[molecule_map_starts_offset];
  super_state_pointers->molecule_map        = &super_statel[molecule_map_offset];
  super_state_pointers->compartment_map     = &super_statel[compartment_map_offset];
  super_state_pointers->molecule_names      = &super_statel[molecule_names_offset];
  super_state_pointers->compartment_names   = &super_statel[compartment_names_offset];
  super_state_pointers->molecules_text       = &super_statec[molecules_text_offset];
  super_state_pointers->compartments_text    = &super_statec[compartments_text_offset];
  return(success);
}
