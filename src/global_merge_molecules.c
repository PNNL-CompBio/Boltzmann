#include "boltzmann_structs.h"
#include "sort_global_molecules.h"
#include "unique_molecules_core.h"
#include "global_merge_molecules.h"
int global_merge_molecules(struct boot_state_struct *boot_state) {
  /*
    Merge local lists of molecules into one list
    Called by: boltzmann_boot
    Calls:     sort_global_molecules, unique_molecules_core
  */
  struct super_state_struct *super_state;
  struct molecule_struct    *molecule_sort_ws;
  struct molecule_struct    *unsorted_molecules;
  struct molecule_struct    *sorted_molecules;
  int64_t align_len;
  int64_t align_mask;
  int64_t global_number_of_molecules;
  int64_t unique_global_molecules;
  int64_t molecule_text_length;
  int64_t num_reaction_files;
  int64_t *molecule_map_starts;
  int64_t *molecule_map;
  char    *molecule_text_ws;
  char    *solvent_string;
  int success;
  int solvent_pos;

  success                    = 1;
  molecule_sort_ws           = boot_state->molecule_sort_ws;
  molecule_text_ws           = boot_state->molecule_text_ws;
  molecule_map_starts        = boot_state->molecule_map_starts;
  molecule_map               = boot_state->molecule_map;
  num_reaction_files         = boot_state->num_reaction_files;
  solvent_string             = boot_state->solvent_string;
  super_state                = boot_state->super_state;
  global_number_of_molecules = super_state->global_number_of_molecules;
  unsorted_molecules = (struct molecule_struct *)&molecule_sort_ws[0];
  sorted_molecules   = (struct molecule_struct *)&molecule_sort_ws[global_number_of_molecules];
  success = sort_global_molecules(&unsorted_molecules,
				  &sorted_molecules,
				  molecule_map_starts,
				  molecule_text_ws,
				  num_reaction_files);
  if (success) {
    /*
      Now we need to extract the list of unique global molecules.    */
    align_len      = boot_state->align_len;
    align_mask     = boot_state->align_mask;
    solvent_string = boot_state->solvent_string;
    success = unique_molecules_core(global_number_of_molecules,
				    sorted_molecules,
				    molecule_text_ws,
				    solvent_string,
				    molecule_map,
				    &unique_global_molecules,
				    &molecule_text_length,
				    &solvent_pos,
				    align_len,align_mask);
  }
  if (success) {
    boot_state->solvent_pos                    = solvent_pos;
    boot_state->sorted_molecules               = sorted_molecules;
    super_state->unique_global_molecules       = unique_global_molecules;
    super_state->molecule_text_length_in_bytes = molecule_text_length;
  }
  return(success);
}
