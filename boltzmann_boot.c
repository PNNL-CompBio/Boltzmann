/* boltzmann_boot.c
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
#include "boot_init.h"
#include "parse_rxn_list_line.h"
#include "boltzmann_init_core.h"
#include "save_and_count_local_state.h"
#include "boot_alloc3.h"
#include "catenate_compartments_and_molecules.h"
#include "global_merge_and_map_compartments.h"
#include "global_merge_molecules.h"
#include "boot_alloc4.h"
#include "condense_strings.h"
#include "fill_meta_data.h"
#include "write_super_state.h"
#include "boltzmann_mmap_superstate.h"
/*
#define DBG_BOLTZMANN_BOOT  
*/
#include "boltzmann_boot.h"
int boltzmann_boot(char *param_file_name, 
		   struct super_state_struct **super_statep) {
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
    compartmtent_text  [compartment_text_length]
----------
    compartment_list_starts [number_of_reaction_files+1]
      These are indexes to the first element per reaction file 
      for the compartments_list vector.
    compartment_list_starts is only needed as workspace, no need to save it.

    compartment_list   [global_number_of_compartments]
    compartment_list is only needed as workspace, no need to save it.

    Called by: boltzmann
    Calls:     alloc0,
               read_params,
	       boot_init,
	       parse_rxn_list_line,
	       boltzmann_init_core,
	       save_and_count_local_state,
	       boot_alloc3,
	       catenate_compartments_and_molecules,
	       global_merge_and_map_compartments,
	       global_merge_molecules,
	       boot_alloc4,
	       condense_strings,
	       fill_meta_data,
	       write_super_state,
	       boltzmann_mmap_superstate,
	       fclose
  */
  struct boot_state_struct  *boot_state;
  struct state_struct *init_state;
  struct state_struct *local_state;
  int64_t mmap_file_len;
  int64_t num_reaction_files;
  int64_t i;
  int64_t *meta_data;
  char    *global_state_filename;
  int success;
  int padi;

  FILE *global_state_fp;
  FILE *lfp;
  /*
    allocate space for the a sample boltzman state struct.
  */
  success = alloc0(&init_state);
  if (success) {
    /*
      Read the input parameters file.
    */
    success = read_params(param_file_name,init_state);
  }
  if (success) {
    /*
      Allocate space for the boot_state struct  and set the boot filenames: 
      rxn_list_file, log file, global state file and work file,
      Also allocate space for the rxn_list_file buffer, size the rxn_list_file, and allocate space
      for the boot meta_data: the superstate_struct, molecule_map_starts, molecule_map, and 
      compartment_list_starts.
    */
    success = boot_init(init_state,&boot_state);
  }
  if (success) {
    num_reaction_files = boot_state->num_reaction_files;
    /*
      Loop over the reaction files.
    */
    for (i=0;i<num_reaction_files;i++) {
      /*
	Parse the reaction list file lines to set the
	reaction_file, input concentrations_file and compartment_sizes_file or sbml_file
	fields of the local_state struct for this reaction group.
      */
      success = parse_rxn_list_line(init_state,boot_state,i);
      if (success) {
	success = boltzmann_init_core(init_state,&local_state);
      }
      if (success) {
	/*
	  Now we need to write the statep block out to a file
	  tracking where we put it for global dictionary determination
	  later.  This stuff needs to be in a routine.
	*/
	success = save_and_count_local_state(local_state,boot_state,i);
      } /* end if we could open and read reaction file */
    } /* end for (i...) */
  }
  if (success) {
    /*
      Need to close and reopen tmp_state_fp to flush the file buffers
      out to the file so that it may be read. and we need to
      allocate space to catenate and sort the compartment and 
      molecule names for all the reactions.
    */
    success = boot_alloc3(boot_state);
  }
  if (success) {
    /*
      Now loop over the local reaction workspaces building
      compartment_text_ws, and molecule_text_ws strings
      as well as the first half of the sort workspaces
      for molecules and compartments.
    */
    success = catenate_compartments_and_molecules(boot_state);
  }
  /*
    at this juncture 
    compartment_sort_ws is a vector of compartment_struct's of length 
    2*global_number_of_compartments, the first half is filled
    with compartments where the c_index field is the global index 
    of the compartment, the g_index field is the reaction file number
    and the string field is an offset into the compartment_text_ws array.
    
    Similarly
    molecule_sort_ws in a vector of molecule_struct's of length
    2*global_number_of_molecules, the first half is filled with
    with molecules where the m_index field is the global index of the
    molecule, the c_index field is the local compartment number
    within the reaction file, the g_index field is the reaction file 
    number, and the string field is an offset into the molecule_text_ws
    array.
  */
  /*
    So now we merge the local compartments, remove duplicates and
    map the local compartmnet numbers to global compartment numbers
    in all of the local molecules structs.
  */
  if (success) {
    success = global_merge_and_map_compartments(boot_state);
  }
  if (success) {
    /*
      So now we sort the molecules, within a reaction file they
      are already sorted so its just a matter of merging these lists.
    */
    success = global_merge_molecules(boot_state);
  }
  /*
    Now we need to allocate space for the
    compartment_names and the molecule_names character pointer arrays,
    and the condensed molecule_text and compartment text arrays,
    and the io_buffer.
  */
  if (success) {
    success = boot_alloc4(boot_state);
  }
  if (success) {
    success = condense_strings(boot_state);
  }
  if (success) {
    success = fill_meta_data(boot_state);
  }
  if (success) {
    success = write_super_state(boot_state); 
  }
  if (success) {
    global_state_fp = boot_state->global_state_fp;
    meta_data       = boot_state->meta_data;
    if (global_state_fp) {
      fclose(global_state_fp);
    }
    mmap_file_len   = meta_data[1];
    global_state_filename = boot_state->global_state_file;
    success = boltzmann_mmap_superstate(global_state_filename,
					mmap_file_len,
					super_statep);
  }
  return(success);
}
