#include "boltzmann_structs.h"
#include "condense_strings.h"
int condense_strings(struct boot_state_struct *boot_state) {
  /*
    Build a condensed list of compartment and molecule names with 
    duplicates removed, and addust the local compartment and molecle
    string pointers to point into the condenesed string arrays.
    Called by: boltzmann_boot
    Calls:     strlen, strcpy
  */
  struct super_state_struct *super_state;
  struct compartment_struct *compartments;
  struct molecule_struct    *molecules;
  int64_t cmpt_size;
  int64_t mlcl_size;
  int64_t pad_size;
  int64_t align_len;
  int64_t align_mask;
  int64_t *compartment_names;
  int64_t *molecule_names;
  int64_t compartment_pos;
  int64_t molecule_pos;
  int64_t unique_global_compartments;
  int64_t unique_global_molecules;
  int64_t i;
  int64_t one_l;
  char    *compartment_text_ws;
  char    *molecule_text_ws;
  char    *cstring;
  char    *mstring;
  char    *compartments_text;
  char    *molecules_text;
  int success;
  int padi;
  success = 1;
  align_len = boot_state->align_len;
  align_mask = boot_state->align_mask;
  compartment_text_ws  = boot_state->compartment_text_ws;
  compartments_text    = boot_state->compartments_text;
  compartment_names    = boot_state->compartment_names;
  compartments         = boot_state->sorted_compartments;
  molecule_text_ws     = boot_state->molecule_text_ws;
  molecules_text       = boot_state->molecules_text;
  molecule_names       = boot_state->molecule_names;
  molecules            = boot_state->sorted_molecules;
  super_state          = boot_state->super_state;
  unique_global_compartments = super_state->unique_global_compartments;
  unique_global_molecules    = super_state->unique_global_molecules;
  one_l             = (int64_t)1;
  compartment_pos = (int64_t)0;
  for (i=0;i<unique_global_compartments;i++) {
    cstring   = (char*)&compartment_text_ws[compartments->string];
    cmpt_size = (int64_t)strlen(cstring) + one_l;
    pad_size = (align_len  - (cmpt_size & align_mask)) & align_mask;
    strcpy((char*)&(compartments_text[compartment_pos]),cstring);
    compartment_names[i] = compartment_pos;
    compartments->string = compartment_pos;
    compartment_pos += (cmpt_size + pad_size);
    compartments++; /* Caution address arithmetic */
  }
  molecule_pos = (int64_t)0;
  for (i=0;i<unique_global_molecules;i++) {
    mstring = (char*)&molecule_text_ws[molecules->string];
    mlcl_size = (int64_t)strlen(mstring) + one_l;
    pad_size = (align_len  - (mlcl_size & align_mask)) & align_mask;
    strcpy((char*)&(molecules_text[molecule_pos]),mstring);
    molecule_names[i] = molecule_pos;
    molecules->string = molecule_pos;
    molecule_pos += (mlcl_size + pad_size);
    molecules++; /* Caution address arithmetic */
  }
  boot_state->compartment_text_length = compartment_pos;
  boot_state->molecule_text_length = molecule_pos;
  return(success);
}
