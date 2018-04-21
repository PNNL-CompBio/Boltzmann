#include "boltzmann_structs.h"

#include "sort_global_compartments.h"
#include "unique_compartments_core.h"

#include "global_merge_and_map_compartments.h"
int global_merge_and_map_compartments(struct boot_state_struct *boot_state) {
  /*
    Sort the compartments globally by merging local compartment lists,
    Remove duplicates,
    Map local compartments to global compartments in all molecule structs.
    Called by: boltzmann_boot
    Calls:     sort_global_compartments, unique_compartments_core,
  */
  struct super_state_struct *super_state;
  struct compartment_struct *compartment_sort_ws;
  struct compartment_struct *unsorted_compartments;
  struct compartment_struct *sorted_compartments;
  struct molecule_struct    *molecule_sort_ws;
  struct molecule_struct    *molecules;
  int64_t num_reaction_files;
  int64_t global_number_of_compartments;
  int64_t global_number_of_molecules;
  int64_t align_len;
  int64_t align_mask;
  int64_t unique_global_compartments;
  int64_t compartment_text_length;
  int64_t i;
  int64_t *compartment_list;
  int64_t *compartment_list_starts;
  int64_t *compartment_map;
  char    *compartment_text_ws;
  
  int success;
  int padi;
  success                 = 1;
  compartment_sort_ws     = boot_state->compartment_sort_ws;
  compartment_text_ws     = boot_state->compartment_text_ws;
  molecule_sort_ws        = boot_state->molecule_sort_ws;
  num_reaction_files      = boot_state->num_reaction_files;
  compartment_list_starts = boot_state->compartment_list_starts;
  compartment_list        = boot_state->compartment_list;
  compartment_map         = boot_state->compartment_map;
  super_state             = boot_state->super_state;
  global_number_of_compartments = super_state->global_number_of_compartments;
  global_number_of_molecules    = super_state->global_number_of_molecules;

  unsorted_compartments = (struct compartment_struct *)&(compartment_sort_ws[0]);
  sorted_compartments   = (struct compartment_struct *)&(compartment_sort_ws[global_number_of_compartments]);
  success = sort_global_compartments(&unsorted_compartments,
				     &sorted_compartments,
				     compartment_list_starts,
				     compartment_text_ws,
				     num_reaction_files);
  if (success) {
    /*
      Now we need to extract the list of unique global compartments.
    */
    align_len      = boot_state->align_len;
    align_mask     = boot_state->align_mask;
    success = unique_compartments_core(global_number_of_compartments,
				       sorted_compartments,
				       compartment_text_ws,
				       compartment_list,
				       &unique_global_compartments,
				       &compartment_text_length,
				       align_len,align_mask);
  }
  if (success) {
    boot_state->sorted_compartments               = sorted_compartments;
    super_state->unique_global_compartments       = unique_global_compartments;
    super_state->compartment_text_length_in_bytes = compartment_text_length;
    /*
      Replace the local compartment indices with global compartment
      indices in the unsorted_molecules.
    */
    molecules = (struct molecule_struct *)&molecule_sort_ws[0];
    for (i=0;i<global_number_of_molecules;i++) {
      /*
	compartment_list is a list of all unique compartment names within
	each reaction file by reaction file.
	molecules->g_index is the reaction_file number.
	compartment_list_starts[molecules->g_index] points to where
	in the compartment_list is to start, and molecules->c_index gives
	the offset from the start.
      */
      molecules->c_index = compartment_list[compartment_list_starts[molecules->g_index] + molecules->c_index];
      compartment_map[i] = molecules->c_index;
      molecules++; /* Caution address arithmetic */
    }
  }
  return(success);
}
