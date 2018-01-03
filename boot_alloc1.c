#include "boltzmann_structs.h"
#include "boot_alloc1.h"
int boot_alloc1(struct boot_state_struct *boot_state) {
  /*
    Allocate space for the fixed length fields of boot_state:
    boot_log_file, global_state_file, boot_work_file, rxn_list_buffer,
    and page_fill.
    Called by: boot_init
    Calls:     calloc,fprintf,fflush
  */
  int64_t ask_for;
  int64_t one_l;
  int64_t page_size;
  int64_t filename_len;
  int64_t rxn_list_buffer_len;
  int64_t *page_fill;
  char *boot_log_file;
  char *rxn_list_file;
  char *global_state_file;
  char *boot_work_file;
  char *rxn_list_buffer;
  char *solvent_string;

  int success;
  int padi;

  success = 1;
  one_l   = (int64_t)1;
  page_size = boot_state->page_size;
  filename_len = boot_state->filename_len;
  rxn_list_buffer_len = boot_state->rxn_list_buffer_len;
  /* 
     Allocate space for the boot log file.
  */
  ask_for = filename_len;
  boot_log_file = (char *)calloc(one_l,ask_for);
  if (boot_log_file == NULL) {
    fprintf(stderr,"boot_alloc1: Error, could not allocate "
	    "%ld bytes for boot_log_file name\n",ask_for);
    fflush(stderr);
    success = 0;
  }
  /*
    Allocate space for the rxn_list_file.
  */
  if (success) {
    boot_state->boot_log_file = boot_log_file;
    rxn_list_file = (char *)calloc(one_l,ask_for);
    if (rxn_list_file == NULL) {
      fprintf(stderr,"boot_alloc1: Error, could not allocate "
	      "%ld bytes for rxn_list_file name\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  /*
    Allocate space for the global_state_file.
  */
  if (success) {
    boot_state->rxn_list_file = rxn_list_file;
    global_state_file = (char *)calloc(one_l,ask_for);
    if (global_state_file == NULL) {
      fprintf(stderr,"boot_alloc1: Error, could not allocate "
	      "%ld bytes for global_state_file name\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  /*
    Allocate space for the boot_work_file.
  */
  if (success) {
    boot_state->global_state_file = global_state_file;
    boot_work_file = (char *)calloc(one_l,ask_for);
    if (boot_work_file == NULL) {
      fprintf(stderr,"boot_alloc1: Error, could not allocate "
	      "%ld bytes for booto_work_file name\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  /*
    Allocate space for the rxn_list_buffer.
  */
  if (success) {
    boot_state->boot_work_file = boot_work_file;
    ask_for = rxn_list_buffer_len;
    rxn_list_buffer = (char*)calloc(one_l,ask_for);
    if (rxn_list_buffer == NULL) {
      fprintf(stderr,"boot_alloc1: Error, could not allocate "
	      "%ld bytes for rxn_list_buffer\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  /*
    Allocate space for the page_fill buffer.
  */
  if (success) {
    boot_state->rxn_list_buffer = rxn_list_buffer;
    ask_for = page_size;
    page_fill = (int64_t*)calloc(one_l,ask_for);
    if (page_fill == NULL) {
      fprintf(stderr,"boot_alloc1: Error, could not allocate "
	      "%ld bytes for page_fill\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  if (success) {
    boot_state->page_fill = page_fill;
    ask_for = (int64_t)64;
    solvent_string = (char*)calloc(one_l,ask_for);
    if (solvent_string == NULL) {
      fprintf(stderr,"boot_alloc1: Error, could not allocate "
	      "%ld bytes for solvent_string\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  if (success) {
    boot_state->solvent_string = solvent_string;
  }
  return(success);
}


