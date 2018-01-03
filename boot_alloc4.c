#include "boltzmann_structs.h"
#include "boot_alloc4.h"
int boot_alloc4 (struct boot_state_struct *boot_state) {
  /*
    Allocate space for the compressed molecule_text,
    and the  number to name translation table for molecules, molecule_names,
    and the  number to name translation table for compartments,
    compartment_names, and the io_buff; 

    Called by: boltzmann_boot
    Calls:     calloc, fprintf, fflush
  */
  struct super_state_struct *super_state;
  char *molecules_text;
  char *compartments_text;
  int64_t *molecule_names;
  int64_t *compartment_names;
  char    *io_buff;
  int64_t one_l;
  int64_t ask_for;
  int64_t int64_t_size;
  int64_t molecule_text_length;
  int64_t compartment_text_length;
  int64_t unique_global_molecules;
  int64_t unique_global_compartments;
  int64_t page_size;
  int64_t io_buff_size_in_pages;

  int success;
  int padi;
  
  FILE *lfp;
  FILE *efp;
  
  success                    = 1;
  one_l                      = (int64_t)1;
  int64_t_size               = (int64_t)sizeof(int64_t);
  lfp                        = boot_state->lfp;
  super_state                = boot_state->super_state;
  molecule_text_length       = super_state->molecule_text_length_in_bytes;
  compartment_text_length    = super_state->compartment_text_length_in_bytes;
  unique_global_molecules    = super_state->unique_global_molecules;
  unique_global_compartments = super_state->unique_global_compartments;
  page_size                  = boot_state->page_size;
  io_buff_size_in_pages      = boot_state->io_buff_size_in_pages;
  
  molecules_text = (char*)calloc(one_l,molecule_text_length);
  if (molecules_text == NULL) {
    success = 0;
    if (lfp) {
      fprintf(lfp,"boot_alloc4: Error could not allocate %ld "
	      "bytes for molecules_text\n",molecule_text_length);
      fflush(lfp);
    }
  }
  if (success) {
    boot_state->molecules_text = molecules_text;
    compartments_text = (char*)calloc(one_l,compartment_text_length);
    if (compartments_text == NULL) {
      success = 0;
      if (lfp) {
    	fprintf(lfp,"boot_alloc4: Error could not allocate %ld "
    		"bytes for compartments_text\n",compartment_text_length);
    	fflush(lfp);
      }
    }
  }
  if (success) {
    boot_state->compartments_text = compartments_text;
    ask_for = int64_t_size * unique_global_molecules;
    molecule_names = (int64_t*)calloc(one_l,ask_for);
    if (molecule_names == NULL) {
      success = 0;
      if (lfp) {
    	fprintf(lfp,"boot_alloc4: Error could not allocate %ld "
    		"bytes for molecule_names\n",ask_for);
    	fflush(lfp);
      }
    }
  }
  if (success) {
    boot_state->molecule_names = molecule_names;
    ask_for = int64_t_size * unique_global_compartments;
    compartment_names = (int64_t*)calloc(one_l,ask_for);
    if (compartment_names == NULL) {
      success = 0;
      if (lfp) {
    	fprintf(lfp,"boot_alloc4: Error could not allocate %ld "
    		"bytes for compartment_names\n",ask_for);
    	fflush(lfp);
      }
    }
  }
  if (success) {
    boot_state->compartment_names = compartment_names;
    ask_for = io_buff_size_in_pages * page_size;
    io_buff = (char *)calloc(one_l,ask_for);
    if (io_buff == NULL) {
      if (lfp) {
	fprintf(lfp,"boot_alloc4: Unable to allocate %ld bytes for io_buff. Will try single page size.\n",ask_for);
	fflush(lfp);
      }
      boot_state->io_buff_size_in_pages = 1;
      ask_for = page_size;
      io_buff = (char *)calloc(one_l,ask_for);
      if (io_buff == NULL) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"boot_alloc4: Unable to allocate %ld bytes for io_buff.\n",ask_for);
	  fflush(lfp);
	}
      }
    }
  }
  if (success) {
    boot_state->io_buff = io_buff;
  }
  return(success);
}
