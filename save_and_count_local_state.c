#include "boltzmann_structs.h"

#include "save_and_count_local_state.h"
int save_and_count_local_state(struct state_struct *local_state,struct boot_state_struct *boot_state, int rxn_num) {
  /*
    Save the local state data to the work file, save the local state length 
    and offset in the work file, and record the molecule text length, 
    number unique molecules, compartment text length, and number of unique 
    compartments, add to the partial sums of unique molecules and unique 
    compartments and update the work file position variable.
    Called by: boltzmann_boot
    Calls:     fprintf, fflush
  */
  int64_t cur_size;
  int64_t fill_size;
  int64_t page_size;
  int64_t page_mask;
  int64_t one_l;
  int64_t nw;
  int64_t offset_size_pair_pos;
  int64_t nunique_molecules;
  int64_t nunique_compartments;
  int64_t *page_fill;
  int64_t *state_offset_size;
  int64_t *molecule_map_starts;
  int64_t *compartment_list_starts;
  char    *rxn_file_name;
  int boot_work_fd;
  int success;
  FILE *lfp;
  FILE *efp;
  success      	    	  = 1;
  one_l        	    	  = (int64_t)1;
  lfp          	    	  = boot_state->lfp;
  boot_work_fd 	    	  = boot_state->boot_work_fd;
  page_fill    	    	  = boot_state->page_fill;
  state_offset_size 	  = boot_state->state_offsets_sizes;
  molecule_map_starts     = boot_state->molecule_map_starts;
  compartment_list_starts = boot_state->compartment_list_starts;
  page_size               = boot_state->page_size;
  page_mask               = page_size - one_l;
  rxn_file_name     	  = local_state->reaction_file;
  boot_state->molecule_text_length    += local_state->molecule_text_length;
  boot_state->compartment_text_length += local_state->compartment_text_length;
  cur_size = local_state->state_length;
  if (cur_size > boot_state->maximum_state_size) {
    boot_state->maximum_state_size = cur_size;
  }
  if ((boot_state->minimum_state_size == 0) || (cur_size < boot_state->minimum_state_size)) {
    boot_state->minimum_state_size = cur_size;
  }
  fill_size = (page_size - (cur_size & page_mask)) & page_mask;
  nw = write(boot_work_fd,local_state,cur_size);
  if (nw != cur_size) {
    success = 0;
    if (lfp) {
      fprintf(lfp,"save_and_count_local_state: Error writing data for %s state\n",
	      rxn_file_name);
      fflush(lfp);
    }
  }
  if (success) {
    if (fill_size > 0) {
      nw = write(boot_work_fd,(char *)page_fill,fill_size);
      if (nw != fill_size) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"boltzmann_boot: Error writing fill data for %s state\n",
		  rxn_file_name);
	  fflush(lfp);
	}
      }
    }
  }
  if (success) {
    offset_size_pair_pos                      = rxn_num + rxn_num;
    state_offset_size[offset_size_pair_pos]   = boot_state->work_offset;
    state_offset_size[offset_size_pair_pos+1] = cur_size;
    nunique_molecules = local_state->nunique_molecules;
    nunique_compartments = local_state->nunique_compartments;
    molecule_map_starts[rxn_num+1] = molecule_map_starts[rxn_num] + nunique_molecules;
    compartment_list_starts[rxn_num+1] = compartment_list_starts[rxn_num] + nunique_compartments;
    if (nunique_molecules > boot_state->unique_molecules_max) {
      boot_state->unique_molecules_max = nunique_molecules;
    }
    if (nunique_compartments > boot_state->unique_compartments_max) {
      boot_state->unique_compartments_max = nunique_compartments;
    }
    boot_state->work_offset += cur_size + fill_size;
  }
  return(success);
}
