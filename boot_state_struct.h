#ifndef _BOOT_STATE_STRUCT_H_
#define _BOOT_STATE_STRUCT_H_ 1
struct boot_state_struct {
  struct super_state_struct *super_state;
  struct molecule_struct *molecule_sort_ws;
  struct molecule_struct *sorted_molecules;
  struct compartment_struct *compartment_sort_ws;
  struct compartment_struct *sorted_compartments;
  /*
    Global constants.
  */
  int64_t page_size;
  int64_t page_mask;
  int64_t align_len;
  int64_t align_mask;
  int64_t filename_len;
  int64_t rxn_list_buffer_len;
  int64_t num_reaction_files;
  int64_t meta_data_size;
  int64_t global_number_of_compartments;
  int64_t global_number_of_molecules;
  int64_t io_buff_size_in_pages;
  /*
    local state sizes.
  */
  int64_t state_offsets_sizes_offset;
  int64_t molecule_map_starts_offset;
  int64_t maximum_state_size;
  int64_t minimum_state_size;
  /*
    Global_position counters.
  */
  int64_t work_offset;
  int64_t molecule_text_length;
  int64_t compartment_text_length;
  int64_t unique_molecules_max;
  int64_t unique_compartments_max;
  /*
    Local pointers
  */
  int64_t *meta_data;
  int64_t *state_offsets_sizes; /* length = 2*number_of_reaction_files */
  int64_t *molecule_map_starts; /* length = number_of_reaction_files+1 */
  int64_t *molecule_map;        /* length = global_number_of_molecules */
  int64_t *compartment_map;     /* length = global_number_of_molecules */
  int64_t *molecule_names;      /* length = unique_global_number_of_molecules */
  int64_t *compartment_names;   /* length = unique_global_number_of_compartments */  
  int64_t *page_fill;
  int64_t *compartment_list_starts;

  int64_t *compartment_list;

  char    *molecules_text;       /* length = molecule_text_length */
  char    *compartments_text;    /* length = compartment_text_length */
  char    *solvent_string;       /* length = 64. */
  
  char    *rxn_list_file;
  char    *global_state_file;
  char    *boot_log_file;
  char    *boot_work_file;
  char    *rxn_list_buffer;

  char    *molecule_text_ws;
  char    *compartment_text_ws;
  char    *io_buff;

  int     log2_page_size;
  int     boot_work_fd;
  int     global_state_fd;
  int     solvent_pos;
  FILE    *lfp;
  FILE    *boot_work_fp;
  FILE    *global_state_fp;
  FILE    *rxn_list_fp;
}
;
#endif
