#include "boltzmann_structs.h"
#include "boot_alloc2.h"

int boot_alloc2(struct boot_state_struct *boot_state) {
  /*
    Allocate the fields of the boot state that we call the meta data and 
    other fields of the boot_state_struct that depend on the number of reaction files:
    compartment_list_starts, molecule_map_starts, state_offsets_sizes.
    Called by: boot_init
    Calls:     fprintf, fflush, sizeof
  */
  struct super_state_struct sss;
  struct super_state_struct *super_state;
  int64_t ask_for;
  int64_t one_l;
  int64_t super_state_struct_size;
  int64_t int64_t_mask;
  int64_t int64_t_size;
  int64_t num_reaction_files;
  int64_t offset_size_pairs_size;
  int64_t molecule_map_starts_size;
  int64_t compartment_list_starts_size;
  int64_t state_offsets_sizes_offset;
  int64_t molecule_map_starts_offset;
  int64_t molecule_map_offset;
  int64_t meta_data_size;
  int64_t *meta_data;
  int64_t *state_offsets_sizes;
  int64_t *molecule_map_starts;
  int64_t *molecule_map;
  int64_t *compartment_list_starts;
  int log2_int64_t_size;
  int success;
  FILE *lfp;
  FILE *efp;
  success           = 1;
  log2_int64_t_size = 3;
  one_l             = (int64_t)1;
  int64_t_size      = (int64_t)8;
  int64_t_mask      = int64_t_size - one_l;
  
  lfp                = boot_state->lfp;
  num_reaction_files = boot_state->num_reaction_files;
  molecule_map_starts_size = num_reaction_files + one_l;
  compartment_list_starts_size = molecule_map_starts_size;

  ask_for = compartment_list_starts_size * int64_t_size;
  compartment_list_starts = (int64_t *)calloc(one_l,ask_for);
  if (compartment_list_starts == NULL) {
    success = 0;
    if (lfp) {
      fprintf(lfp,"boot_alloc2: Error unable to allocate %lld bytes for compartment_list_starts\n",
	      ask_for);
      fflush(lfp);
    }
  }
  if (success) {
    boot_state->compartment_list_starts = compartment_list_starts;
    super_state_struct_size = ((int64_t)sizeof(sss) + int64_t_mask) >> log2_int64_t_size;
    offset_size_pairs_size = num_reaction_files + num_reaction_files;
    meta_data_size = super_state_struct_size + offset_size_pairs_size + molecule_map_starts_size;
    boot_state->meta_data_size = meta_data_size;
    ask_for = meta_data_size * int64_t_size;
    meta_data = (int64_t *)calloc(one_l,ask_for);
    if (meta_data == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"boot_alloc2: Error unable to allocate %lld bytes for meta_data\n",
		ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    boot_state->meta_data = meta_data;
    boot_state->super_state = (struct super_state_struct*)meta_data;
    state_offsets_sizes_offset = super_state_struct_size;
    molecule_map_starts_offset = state_offsets_sizes_offset + offset_size_pairs_size;
    molecule_map_offset        = molecule_map_starts_offset + molecule_map_starts_size;
    
    state_offsets_sizes   = (int64_t *)&meta_data[state_offsets_sizes_offset];
    molecule_map_starts   = (int64_t *)&meta_data[molecule_map_starts_offset];
    molecule_map          = (int64_t *)&meta_data[molecule_map_offset];
    boot_state->state_offsets_sizes = state_offsets_sizes;
    boot_state->molecule_map_starts = molecule_map_starts;
    boot_state->molecule_map        = molecule_map;
    molecule_map_starts[0]          = (int64_t)0;
    compartment_list_starts[0]      = (int64_t)0;
  }
  return(success);
}
