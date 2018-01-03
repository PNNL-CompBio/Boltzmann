#include "boltzmann_structs.h"
#include "fill_meta_data.h"
int fill_meta_data(struct boot_state_struct *boot_state) {
  /*
    Fill the super_state meta data vector and adjust the offset_size
    pairs for local states.
    Called by: boltzmann_boot
    Calls:
  */
  struct super_state_struct *super_state;
  int64_t *meta_data;
  int64_t *state_offset_size_pair;
  int64_t meta_data_size;
  int64_t molecule_map_offset;
  int64_t compartment_map_offset;
  int64_t molecule_names_offset;
  int64_t molecules_text_offset;
  int64_t compartment_names_offset;
  int64_t compartments_text_offset;
  int64_t global_number_of_molecules;
  int64_t unique_global_molecules;
  int64_t unique_global_compartments;
  int64_t molecule_text_length;
  int64_t compartment_text_length;
  int64_t fill_size;
  int64_t page_size;
  int64_t page_mask;
  int64_t total_length;
  int64_t i;
  int64_t num_reaction_files;
  
  int success;
  int log2_int64_t_size;
  int log2_page_size;

  success             = 1;
  log2_int64_t_size   = 3;
  meta_data_size             = boot_state->meta_data_size;
  meta_data                  = boot_state->meta_data;
  super_state                = boot_state->super_state;
  global_number_of_molecules = super_state->global_number_of_molecules;
  unique_global_molecules    = super_state->unique_global_molecules;
  unique_global_molecules    = super_state->unique_global_compartments;
  compartment_text_length    = boot_state->compartment_text_length;
  molecule_text_length       = boot_state->molecule_text_length;
  page_size                  = boot_state->page_size;
  page_mask                  = boot_state->page_mask;
  log2_page_size             = boot_state->log2_page_size;
  num_reaction_files         = boot_state->num_reaction_files;

  molecule_map_offset        = meta_data_size;
  compartment_map_offset     = meta_data_size + global_number_of_molecules;
  molecule_names_offset      = compartment_map_offset + global_number_of_molecules;
  compartment_names_offset   = molecule_names_offset + unique_global_molecules;
  molecules_text_offset      = (compartment_names_offset + unique_global_compartments)<<log2_int64_t_size;
  compartments_text_offset   = molecules_text_offset + molecule_text_length;
  meta_data_size             = (compartments_text_offset + compartment_text_length);
  fill_size = (page_size - (meta_data_size & page_mask)) & page_mask;
  meta_data_size += fill_size;
  boot_state->meta_data_size = meta_data_size;
  /*
    Fill the meta data array.
  */
  total_length = meta_data_size + boot_state->work_offset;
  super_state->number_of_reaction_files = boot_state->num_reaction_files;
  super_state->total_length_in_bytes = total_length;
  super_state->page_size_in_bytes = page_size;
  super_state->number_of_pages = (total_length >> log2_page_size) +
      ((page_size - (total_length & page_mask)) & page_mask);
  super_state->string_align_len = boot_state->align_len;
  super_state->string_align_mask = boot_state->align_mask;
  super_state->molecule_text_length_in_bytes = molecule_text_length;
  super_state->compartment_text_length_in_bytes = compartment_text_length;
  super_state->state_offsets_sizes_offset_in_words = boot_state->state_offsets_sizes_offset;
  super_state->molecule_map_starts_offset_in_words = boot_state->molecule_map_starts_offset;
  super_state->molecule_map_offset_in_words      = molecule_map_offset;
  super_state->molecule_names_offset_in_words    = molecule_names_offset;
  super_state->compartment_map_offset_in_words   = compartment_map_offset;
  super_state->compartment_names_offset_in_words = compartment_names_offset;
  super_state->molecules_text_offset_in_bytes    = molecules_text_offset;
  super_state->compartments_text_offset_in_bytes = compartments_text_offset;
  super_state->maximum_state_size_in_bytes       = boot_state->maximum_state_size;
  super_state->minimum_state_size_in_bytes       = boot_state->minimum_state_size;
  /*
    Increment the state_offset_size offsets by meta_data_size;
  */
  state_offset_size_pair  = boot_state->state_offsets_sizes;
  for (i=(int64_t)0;i<num_reaction_files;i++) {
    *state_offset_size_pair += meta_data_size;
    state_offset_size_pair += 2; /* Caution address arithmetic */
  }
  return(success);
}
