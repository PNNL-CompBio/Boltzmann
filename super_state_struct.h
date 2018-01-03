/* super_state_struct.h 
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
#ifndef __SUPER_STATE_STRUCT__
#define __SUPER_STATE_STRUCT__ 1

  /*
    Initialize the reactions and data structures for boltzmann.
    This is a variation on boltzmann_init in that it processes
    multiple input reaction files, and builds a single output file of the
    form 
    meta_data[0] 8-byte number of reaction files.
    meta_data[1] 8-byte total  length in bytes.
    meta_data[2] 8-byte page_size in bytes, default val is 4K = 4096
    meta_data[3] 8-byte number of pages = (total_length >> log2(page_size)) + 
                 (page_size - (total_length & (page_size-1))) & (page_size-1))
    meta_data[4] 8-byte align_len for strings, default is 16.
    meta_data[5] 8-byte align_mask             default is 15.
    meta_data[6] 8-byte global number of molecules.
    meta_data[7] 8-byte unique_global_molecules
    meta_data[8] 8-byte molecule_text_length (in bytes)
    meta_data[9] 8-byte global number of compartments. 
    meta_data[10] 8-byte unique_global_compartments 
    meta_data[11] 8-byte compartment_text_length (in bytes)
    meta_data[12] 8 byte state_offsets_sizes_offset (in 8byte words)
    meta_data[13] 8 byte molecule_map_starts_offset; (in 8byte words)
    meta_data[14] 8-byte molecule_map_offset (in 8byte words)
    meta_data[15] 8-byte molecule_names_offset (in 8byte words)
    meta_data[16] 8-byte compartment_map_offset (in 8byte words)
    meta_data[17] 8-byte compartment_names_offset (in 8byte words)
    meta_data[18] 8-byte molecules_text_offset (in bytes)
    meta_data[19] 8-byte compartments_text_offset (in bytes)

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
----------
    compartment_list_starts [number_of_reaction_files+1]
      These are indexes to the first element per reaction file 
      for the compartments_list vector.
    compartment_list_starts is only needed as workspace, no need to save it.

    compartment_list   [global_number_of_compartments]
    compartment_list is only needed as workspace, no need to save it.
*/
struct super_state_struct {
  int64_t number_of_reaction_files;
  int64_t total_length_in_bytes; /* Includes lengths of arrays pointed to. */
  int64_t page_size_in_bytes; /* default val is 4K = 4096 */
  int64_t number_of_pages; 
  /* number_of_pages = (total_length >> log2(page_size)) + 
	(page_size - (total_length & (page_size-1))) & (page_size-1)) 
  */
  int64_t string_align_len; /* default is 16.*/
  int64_t string_align_mask; /* default is 15.*/
  int64_t global_number_of_molecules;
  int64_t unique_global_molecules;
  int64_t molecule_text_length_in_bytes;
  int64_t global_number_of_compartments;
  int64_t unique_global_compartments; 
  int64_t compartment_text_length_in_bytes;
  int64_t state_offsets_sizes_offset_in_words;
  int64_t molecule_map_starts_offset_in_words;
  int64_t molecule_map_offset_in_words;
  int64_t molecule_names_offset_in_words;
  int64_t compartment_map_offset_in_words;
  int64_t compartment_names_offset_in_words;
  int64_t molecules_text_offset_in_bytes;
  int64_t compartments_text_offset_in_bytes;
  int64_t minimum_state_size_in_bytes;
  int64_t maximum_state_size_in_bytes;
}
;
#endif
