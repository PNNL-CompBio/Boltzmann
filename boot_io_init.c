#include "boltzmann_structs.h"
#include "boot_io_init.h"
int boot_io_init(struct state_struct *state, 
		 struct boot_state_struct *boot_state) {
  /*
    Create filenames for and open the boot_log_file, global_state_file and 
    boot_work_file.
    Called by: boot_init
    Calls:     strlen,strncpy,strcpy,fprintf,fflush,fileno
  */
  int success;
  int k;
  int j;
  int i;
  int max_base;
  int padi;
  
  char *rxn_list_file;
  char *boot_log_file;
  char *global_state_file;
  char *boot_work_file;
  char *solvent_string;

  FILE *lfp;
  FILE *boot_work_fp;
  FILE *global_state_fp;
  FILE *efp;

  max_base          = (int)boot_state->filename_len - 7;
  boot_log_file     = boot_state->boot_log_file;
  rxn_list_file     = boot_state->rxn_list_file;
  global_state_file = boot_state->global_state_file;
  boot_work_file    = boot_state->boot_work_file;
  solvent_string    = boot_state->solvent_string;
  strcpy(solvent_string,state->solvent_string);
  strcpy(rxn_list_file,state->reaction_file);
  k             = strlen(rxn_list_file);
  i             = k;
  for (j=k-1;j>=0;j--) {
    if (rxn_list_file[j] == '.') {
      i = j;
      break;
    }
  }
  if (i > max_base) {
    i = max_base;
  }
  strncpy(boot_log_file,rxn_list_file,i);
  strcpy((char*)&boot_log_file[i],".log");
  lfp = fopen(boot_log_file,"w");
  if (lfp == NULL) {
    success = 0;
    fprintf(stderr,"boot_io_init: Unable to open log file, %s, for writing.\n",
	    boot_log_file);
    fflush(stderr);
  }
  if (success) {
    boot_state->lfp = lfp;
    strncpy(boot_work_file,rxn_list_file,i);
    strcpy((char*)&boot_work_file[i],".work");
    boot_work_fp = fopen(boot_work_file,"w");
    if (boot_work_fp == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"boot_io_init: Unable to open work file, %s, for writing.\n",
		boot_work_file);
	fflush(lfp);
      }
    }
  }
  if (success) {
    boot_state->boot_work_fp = boot_work_fp;
    boot_state->boot_work_fd = fileno(boot_work_fp);
    strncpy(global_state_file,rxn_list_file,i);
    strcpy((char*)&global_state_file[i],".state");
    global_state_fp = fopen(global_state_file,"w");
    if (global_state_fp == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"boot_io_init: Unable to open global_state file, %s, for writing.\n",
		global_state_file);
	fflush(lfp);
      }
    }
  }
  if (success) {
    boot_state->global_state_fp = global_state_fp;
    boot_state->global_state_fd = fileno(global_state_fp);
  }
  return(success);
}
