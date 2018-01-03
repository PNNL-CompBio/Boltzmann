/* super_state_pointers_struct.h 
*******************************************************************************
Boltzmann

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
#ifndef __SUPER_STATE_POINTERS_STRUCT__
#define __SUPER_STATE_POINTERS_STRUCT__ 1
struct super_state_pointers_struct {
  /* 
    state_offsets_sizes[number_of_reaction_files]
       offset in bytes, length pairs of the reaction states.

    molecule_map_starts [number_of_reaction_files+1]
      These are indexes to the first element per reaction file 
      for the molecles_map vector.

    molecule_map           [global_number_of_molecules]
    compartment_map        [global_number_of_molecules]

    molecule_names     [unqiue_global_number_of_molecules]
       These will be character offsets relative to the
       beginning of the molecules_text address.
    
    compartment_names  [unique_global_number_of_compartments]
       These will be character offsets relative to the
       beginning of the compartments_text address.

    molecule_text      [molecule_text_length]
    comparmtent_text  [compartment_text_length]
  */
  int64_t *state_offsets_sizes; /* length = 2*number_of_reaction_files */
  int64_t *molecule_map_starts; /* length = number_of_reaction_files+1 */
  int64_t *molecule_map;        /* length = global_number_of_molecules */
  int64_t *compartment_map;     /* length = global_number_of_molecules */
  int64_t *molecule_names;      /* length = unique_global_number_of_molecules */
  int64_t *compartment_names;   /* length = unique_global_number_of_compartments */  
  char    *molecules_text;       /* length = molecule_text_length */
  char    *compartments_text;    /* length = compartment_text_length */
}
;
#endif  
