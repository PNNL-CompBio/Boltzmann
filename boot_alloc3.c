#include "boltzmann_structs.h"

#include "boot_alloc3.h"
int boot_alloc3(struct boot_state_struct *boot_state) {
  /*
    Close and reopen tmp_state_fp to flush the file buffers
    out to the file so that it may be read. and we need to
    Allocate space to catenate and sort the compartment and 
    molecule names for all the reactions.
    Called by: boltzmann_boot
    Calls:     fclose, fopen, fileno, fprintf, fflush, calloc
  */
  struct super_state_struct *super_state;
  struct molecule_struct molecule_instance;
  struct molecule_struct *molecule_sort_ws;
  struct compartment_struct compartment_instance;
  struct compartment_struct *compartment_sort_ws;
  int64_t *compartment_list_starts;
  int64_t *compartment_map_starts;
  int64_t *molecule_map_starts;
  int64_t *compartment_list;
  int64_t *compartment_map;
  int64_t *molecule_map;
  int64_t ask_for;
  int64_t one_l;
  int64_t num_reaction_files;
  int64_t global_number_of_compartments;
  int64_t global_number_of_molecules;
  int64_t molecule_text_length;
  int64_t compartment_text_length;
  int64_t int64_t_size;
  char    *boot_work_file;
  char    *molecule_text_ws;
  char    *compartment_text_ws;
  int success;
  int boot_work_fd;
  
  FILE *boot_work_fp;
  FILE *lfp;
  success                  = 1;
  one_l                    = (int64_t)1;
  int64_t_size             = (int64_t)sizeof(int64_t);
  lfp                      = boot_state->lfp;
  super_state              = boot_state->super_state;
  compartment_list_starts  = boot_state->compartment_list_starts;
  molecule_map_starts      = boot_state->molecule_map_starts;
  num_reaction_files       = boot_state->num_reaction_files;
  boot_work_file           = boot_state->boot_work_file;
  /*
    The following two variables are set in the save_and_count_local_state
    routine.
  */
  molecule_text_length     = boot_state->molecule_text_length;
  compartment_text_length  = boot_state->compartment_text_length;

  /*
    Close and reopen the boot work file for reading.
  */
  boot_work_fp = boot_state->boot_work_fp;
  fclose(boot_work_fp);
  boot_work_fp = fopen(boot_work_file,"r");
  if (boot_work_fp == NULL) {
    if (lfp) {
      fprintf(lfp,"boot_alloc3: Error, could not open "
	    "boot_work_file, %s, for reading\n",boot_work_file);
      fflush(lfp);
      success = 0;
    }
  }
  if (success) {
    boot_work_fd = fileno(boot_work_fp);
    boot_state->boot_work_fp = boot_work_fp;
    boot_state->boot_work_fd = boot_work_fd;
    global_number_of_compartments = compartment_list_starts[num_reaction_files];
    global_number_of_molecules = molecule_map_starts[num_reaction_files];
    super_state->global_number_of_compartments = global_number_of_compartments;
    super_state->global_number_of_molecules    = global_number_of_molecules;
    /*
      Now we need to allocate workspace to store the 
      molecules_text and compartment_text pieces as well as
      double the number of quadruples and triples for each. for
      sorting.
    */
    molecule_text_ws = (char*)calloc(one_l,molecule_text_length);
    if (molecule_text_ws == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"boot_alloc3: Error could not allocate %lld "
    		"bytes for molecule_text_ws\n",molecule_text_length);
    	fflush(lfp);
      }
    }
  }
  if (success) {
    boot_state->molecule_text_ws = molecule_text_ws;
    compartment_text_ws = (char*)calloc(one_l,compartment_text_length);
    if (compartment_text_ws == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"boot_alloc3: Error could not allocate %lld "
	      "bytes for compartment_text_ws\n",compartment_text_length);
	fflush(lfp);
      }
    }
  }
  if (success) {
    boot_state->compartment_text_ws = compartment_text_ws;
    ask_for = (int64_t)(2*sizeof(molecule_instance)) * global_number_of_molecules;
    molecule_sort_ws = (struct molecule_struct *)calloc(one_l,ask_for);
    if (molecule_sort_ws == NULL) {
      success  = 0;
      if (lfp) {
	fprintf(lfp,"boot_alloc3: Error could not allocate %lld "
		"bytes for molecule_sort_ws\n",ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    boot_state->molecule_sort_ws = molecule_sort_ws;
    ask_for = (int64_t)(2*sizeof(compartment_instance)) * 
      global_number_of_compartments;
    compartment_sort_ws = (struct compartment_struct*)calloc(one_l,ask_for);
    if (compartment_sort_ws == NULL) {
      if (lfp) {
	fprintf(lfp,"boot_alloc3: Error could not allocate %lld "
	      "bytes for compartment_sort_ws\n",ask_for);
	fflush(lfp);
	success  = 0;
      }
    }
  }
  if (success) {
    boot_state->compartment_sort_ws = compartment_sort_ws;
    ask_for = (int64_t)(int64_t_size*global_number_of_compartments);
    compartment_list = (int64_t*)calloc(one_l,ask_for);
    if (compartment_list == NULL) {
      success  = 0;
      if (lfp) {
	fprintf(lfp,"boot_alloc3: Error could not allocate %lld "
		"bytes for compartment_list\n",ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    boot_state->compartment_list = compartment_list;
    ask_for = (int64_t)(int64_t_size*global_number_of_molecules);
    compartment_map = (int64_t*)calloc(one_l,ask_for);
    if (compartment_map == NULL) {
      success  = 0;
      if (lfp) {
	fprintf(lfp,"boot_alloc3: Error could not allocate %lld "
	      "bytes for compartment_map\n",ask_for);
	fflush(lfp);
      }
    }
  }	    
  if (success) {
    boot_state->compartment_map = compartment_map;
    ask_for = (int64_t)(int64_t_size*global_number_of_molecules);
    molecule_map = (int64_t*)calloc(one_l,ask_for);
    if (molecule_map == NULL) {
      success  = 0;
      if (lfp) {
	fprintf(lfp,"boot_alloc3: Error could not allocate %lld "
	      "bytes for molecule_map\n",ask_for);
	fflush(lfp);
      }
    }
  }	    
  if (success) {
    boot_state->molecule_map = molecule_map;
  }
  return(success);
}
