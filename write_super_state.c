#include "boltzmann_structs.h"

#include "copy_local_states.h"

#include "write_super_state.h"
int write_super_state(struct boot_state_struct *boot_state) {
  /*
    Write out the self describing global(super)_state file.
    Called by: boltzmann_boot;
    Calls:     copy_local_state,write,fprintf,fflush
  */
  struct super_state_struct *super_state;
  int64_t ntw;
  int64_t nw;
  int64_t molecule_map_offset;
  int64_t global_number_of_molecules;
  int64_t unique_global_molecules;
  int64_t unique_global_compartments;
  int64_t molecule_text_length;
  int64_t compartment_text_length;
  int64_t compartments_text_offset;
  int64_t fill_size;
  int64_t int64_t_size;
  int64_t page_size;
  int64_t page_mask;
  int64_t meta_data_size;
  int64_t *meta_data;
  int64_t *molecule_map;
  int64_t *molecule_names;
  int64_t *compartment_map;
  int64_t *compartment_names;
  int64_t *page_fill;
  char    *molecules_text;
  char    *compartments_text;
  int log2_int64_t_size;
  int success;
  int global_state_fd;
  
  FILE *lfp;
  FILE *efp;

  success             = 1;  
  log2_int64_t_size   = 3;
  int64_t_size        = (int64_t)8;
  lfp                 = boot_state->lfp;
  super_state         = boot_state->super_state;
  global_state_fd     = boot_state->global_state_fd;
  meta_data           = boot_state->meta_data;
  page_size           = boot_state->page_size;
  page_mask           = boot_state->page_mask;
  molecule_map_offset = super_state->molecule_map_offset_in_words;
  global_number_of_molecules    = boot_state->global_number_of_molecules;
  molecule_map                  = boot_state->molecule_map;
  unique_global_molecules       = super_state->unique_global_molecules;
  molecule_names                = boot_state->molecule_names;
  compartment_map               = boot_state->compartment_map;
  unique_global_compartments    = super_state->unique_global_compartments;
  compartment_names             = boot_state->compartment_names;
  molecule_text_length          = boot_state->molecule_text_length;
  molecules_text                = boot_state->molecules_text;
  compartment_text_length       = boot_state->compartment_text_length;
  compartments_text             = boot_state->compartments_text;
  page_fill                     = boot_state->page_fill;
  compartments_text_offset      = super_state->compartments_text_offset_in_bytes;

  ntw = molecule_map_offset << log2_int64_t_size;
  nw = write(global_state_fd,(void*)&meta_data[0],ntw);
  if (nw != ntw) {
    success = 0;
    if (lfp) {
      fprintf(lfp,"write_super_state: Error writing meta data to global state file\n");
      fflush(lfp);
    }
  }
  if (success) {
    ntw = global_number_of_molecules << log2_int64_t_size;
    nw = write(global_state_fd,(void*)&molecule_map[0],ntw);
    if (nw!=ntw) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"write_super_state: Error writing molecule_map to global state file\n");
	fflush(lfp);
      }
    }
  }
  if (success) {
    ntw = global_number_of_molecules << log2_int64_t_size;
    nw = write(global_state_fd,(void*)&compartment_map[0],ntw);
    if (nw!=ntw) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"write_super_state: Error writing compartment_map to global state file\n");
	fflush(lfp);
      }
    }
  }  
  if (success) {
    ntw = unique_global_molecules << log2_int64_t_size;
    nw = write(global_state_fd,(void*)&molecule_names[0],ntw);
    if (nw!=ntw) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"write_super_state: Error writing molecule_names to global state file\n");
	fflush(lfp);
      }
    }
  }  
  if (success) {
    ntw = unique_global_compartments << log2_int64_t_size;
    nw = write(global_state_fd,(int64_t*)&compartment_names[0],ntw);
    if (nw!=ntw) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"write_super_state: Error writing compartment_names to global state file\n");
	fflush(lfp);
      }
    }
  }  
  if (success) {
    ntw = molecule_text_length;
    nw = write(global_state_fd,(void*)molecules_text,ntw);
    if (nw!=ntw) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"write_super_state: Error writing molecule_text to global state file\n");
	fflush(lfp);
      }
    }
  }  
  if (success) {
    ntw = compartment_text_length;;
    nw = write(global_state_fd,(void*)compartments_text,ntw);
    if (nw!=ntw) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"write_super_state: Error writing compartment_text to global state file\n");
	fflush(lfp);
      }
    }
  }  
  if (success) {
    meta_data_size = compartments_text_offset + compartment_text_length;
    fill_size = (page_size - (meta_data_size & page_mask)) & page_mask;
    ntw = fill_size;
    nw = write(global_state_fd,(void*)&page_fill[0],ntw);
    if (nw!=ntw) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"write_super_state: Error writing pagefill to global state file\n");
	fflush(lfp);
      }
    }
  }  
  if (success) {
    success = copy_local_states(boot_state);
  }
  return(success);
}
