#include "boltzmann_structs.h"
#include "create_output_filenames.h"
int create_output_filenames(struct state_struct *state) {
  /*
    Create output file names for the counts_out, rxn_likelihood,
    free_energey, boundary_flux, concs_in, cmpts_dat rxns_dat, id_name,
    restart, rxn_view and rxn_mat files. In the process set the
    rxn_filename_base_length field of state.
    Called by: io_size_init
    Calls:     strlen,strncpy,strcpy,fprintf,fflush
  */
  char *rxn_filename;
  char *init_conc_filename;
  char *log_filename;
  char *output_filename;
  char *counts_out_filename;
  char *ode_concs_filename;
  char *rxn_lklhd_filename;
  char *free_energy_filename;
  char *restart_filename;
  char *rxn_view_filename;
  char *bndry_flux_filename;
  char *rxn_echo_filename;
  char *rxn_mat_filename;
  char *dg0ke_filename;
  char *dictionary_filename;
  int64_t rxn_filename_base_length;
  int rxn_filename_length;
  int success;
  int i;
  int no_slash_or_dot;

  success              = 1;
  rxn_filename         = state->reaction_file;
  init_conc_filename   = state->init_conc_file;
  log_filename         = state->log_file;
  output_filename      = state->output_file;
  counts_out_filename  = state->counts_out_file;
  ode_concs_filename   = state->ode_concs_file;
  rxn_lklhd_filename   = state->rxn_lklhd_file;
  free_energy_filename = state->free_energy_file;
  restart_filename     = state->restart_file;
  rxn_view_filename    = state->rxn_view_file;
  bndry_flux_filename  = state->bndry_flux_file;
  rxn_echo_filename    = state->rxn_echo_file;
  rxn_mat_filename     = state->rxn_mat_file;
  dg0ke_filename       = state->dg0ke_file;
  dictionary_filename  = state->dictionary_file;

  rxn_filename_length = strlen(rxn_filename);
  no_slash_or_dot = 1;
  rxn_filename_base_length = -1;
  for (i=rxn_filename_length-1;((i>=1) && no_slash_or_dot);i--) {
    if (rxn_filename[i] == '.') {
      rxn_filename_base_length = i;
      no_slash_or_dot = 0;
    } else {
      if (rxn_filename[i] == '/') {
	no_slash_or_dot = 0;
      }
    }
  }
  if (rxn_filename_base_length < 1) {
    rxn_filename_base_length = rxn_filename_length;
    /*
      check for a reaction file name of "." or "./",
      if so these will flag as an error.
    */
    if (rxn_filename_length < 3) {
      if ((strcmp(rxn_filename,"./") == 0)  || (strcmp(rxn_filename,".") == 0) || (rxn_filename_length < 1)) {
	success = 0;
	fprintf(stderr,"create_output_filenames: Error invalid rxn_file was %s\n",rxn_filename);
	fflush(stderr);
      }
    }
  }
  if (success) {
    if ((rxn_filename_base_length + 11) > state->max_filename_len) {
      rxn_filename_base_length = state->max_filename_len - 12;
      /*
	Probably should print a warning mewssage here.
      */
      fprintf(stderr,"create_output_filenames: Warning, truncated rxn_filename_base to %ld characters\n",rxn_filename_base_length);
      fflush(stderr);
    }
    state->rxn_filename_base_length = rxn_filename_base_length;
    if (init_conc_filename[0] == '\0') {
      strncpy(init_conc_filename,rxn_filename,rxn_filename_base_length);
      strcpy((char*)&init_conc_filename[rxn_filename_base_length],".concs");
    }
    if (log_filename[0] == '\0') {
      strncpy(log_filename,rxn_filename,rxn_filename_base_length);
      strcpy((char*)&log_filename[rxn_filename_base_length],".log");
    }
    if (output_filename[0] == '\0') {
      strncpy(output_filename,rxn_filename,rxn_filename_base_length);
      strcpy((char*)&output_filename[rxn_filename_base_length],".out");
    }
    if (counts_out_filename[0] == '\0') {
      strncpy(counts_out_filename,rxn_filename,rxn_filename_base_length);
      strcpy((char*)&counts_out_filename[rxn_filename_base_length],".count");
    }
    if (ode_concs_filename[0] == '\0') {
      strncpy(ode_concs_filename,rxn_filename,rxn_filename_base_length);      
      strcpy((char*)&ode_concs_filename[rxn_filename_base_length],".ode_concs");
    }
    if (rxn_lklhd_filename[0] == '\0') {
      strncpy(rxn_lklhd_filename,rxn_filename,rxn_filename_base_length);
      strcpy((char*)&rxn_lklhd_filename[rxn_filename_base_length],".lklhd");
    }
    if (free_energy_filename[0] == '\0') {
      strncpy(free_energy_filename,rxn_filename,rxn_filename_base_length);
      strcpy((char*)&free_energy_filename[rxn_filename_base_length],".fe");
    }
    if (restart_filename[0] == '\0') {
      strncpy(restart_filename,rxn_filename,rxn_filename_base_length);
      strcpy((char*)&restart_filename[rxn_filename_base_length],".rstrt");
    }
    if (rxn_view_filename[0] == '\0') {
      strncpy(rxn_view_filename,rxn_filename,rxn_filename_base_length);
      strcpy((char*)&rxn_view_filename[rxn_filename_base_length],".view");
    }
    if (bndry_flux_filename[0] == '\0') {
      strncpy(bndry_flux_filename,rxn_filename,rxn_filename_base_length);
      strcpy((char*)&bndry_flux_filename[rxn_filename_base_length],".flux");
    }
    if (rxn_echo_filename[0] == '\0') {
      strncpy(rxn_echo_filename,rxn_filename,rxn_filename_base_length);
      strcpy((char*)&rxn_echo_filename[rxn_filename_base_length],".echo");
    }
    if (rxn_mat_filename[0] == '\0') {
      strncpy(rxn_mat_filename,rxn_filename,rxn_filename_base_length);
      strcpy((char*)&rxn_mat_filename[rxn_filename_base_length],".mat");
    }
    if (dg0ke_filename[0] == '\0') {
      strncpy(dg0ke_filename,rxn_filename,rxn_filename_base_length);
      strcpy((char*)&dg0ke_filename[rxn_filename_base_length],".dg0ke");
    }
    if (dictionary_filename[0] == '\0') {
      strncpy(dictionary_filename,rxn_filename,rxn_filename_base_length);
      strcpy((char*)&dictionary_filename[rxn_filename_base_length],".dict");
    }
  }
  return(success);
}
