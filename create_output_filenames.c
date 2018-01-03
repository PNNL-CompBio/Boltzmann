#include "boltzmann_structs.h"
#include "create_output_filenames.h"
int create_output_filenames(struct state_struct *state) {
  /*
    Create output file names for the counts_out, rxn_likelihood,
    free_energey, boundary_flux, concs_in, cmpts_dat rxns_dat, id_name,
    restart, rxn_view and rxn_mat files. In the process set the
    output_filename_base_length field of state.
    Called by: io_size_init
    Calls:     strlen,strncpy,strcpy,fprintf,fflush
  */
  char *rxn_filename;
  char *init_conc_filename;
  char *log_filename;
  char *output_filename;
  char *counts_out_filename;
  char *concs_out_filename;
  char *ode_concs_filename;
  char *ode_dconcs_filename;
  char *ode_lklhd_filename;
  char *ode_bflux_filename;
  char *net_lklhd_filename;
  char *nl_bndry_flx_filename;
  char *rxn_lklhd_filename;
  char *free_energy_filename;
  char *restart_filename;
  char *rxn_view_filename;
  char *bndry_flux_filename;
  char *rxn_echo_filename;
  char *rxn_mat_filename;
  char *dg0ke_filename;
  char *dictionary_filename;
  char *aux_data_filename;
  char *ode_counts_filename;
  char *ode_sens_filename;
  char *ode_dsens_filename;
  int64_t output_filename_base_length;
  int output_filename_length;
  int success;
  int i;
  int no_slash_or_dot;
  int output_filename_empty;
  int padi;
  success              = 1;
  rxn_filename         = state->reaction_file;
  init_conc_filename   = state->init_conc_file;
  log_filename         = state->log_file;
  output_filename      = state->output_file;
  counts_out_filename  = state->counts_out_file;
  concs_out_filename   = state->concs_out_file;
  ode_concs_filename   = state->ode_concs_file;
  ode_counts_filename  = state->ode_counts_file;
  net_lklhd_filename   = state->net_lklhd_file;
  nl_bndry_flx_filename = state->nl_bndry_flx_file;
  rxn_lklhd_filename   = state->rxn_lklhd_file;
  free_energy_filename = state->free_energy_file;
  restart_filename     = state->restart_file;
  rxn_view_filename    = state->rxn_view_file;
  bndry_flux_filename  = state->bndry_flux_file;
  rxn_echo_filename    = state->rxn_echo_file;
  rxn_mat_filename     = state->rxn_mat_file;
  dg0ke_filename       = state->dg0ke_file;
  ode_dconcs_filename  = state->ode_dconcs_file;
  ode_lklhd_filename   = state->ode_lklhd_file;
  ode_bflux_filename   = state->ode_bflux_file;
  aux_data_filename    = state->aux_data_file;
  ode_sens_filename    = state->ode_sens_file;
  ode_dsens_filename   = state->ode_dsens_file;

  dictionary_filename  = state->dictionary_file;


  output_filename_empty = 0;
  output_filename_length = strlen(output_filename);
  if (output_filename_length == 0) {
    output_filename_empty = 1;
    /*
      No output file was specified so use input file name.
    */
    strcpy(output_filename,rxn_filename);
  }
  no_slash_or_dot = 1;
  output_filename_base_length = -1;
  for (i=output_filename_length-1;((i>=1) && no_slash_or_dot);i--) {
    if (output_filename[i] == '.') {
      output_filename_base_length = i;
      no_slash_or_dot = 0;
    } else {
      if (output_filename[i] == '/') {
	no_slash_or_dot = 0;
      }
    }
  }
  if (output_filename_base_length < 1) {
    output_filename_base_length = output_filename_length;
    /*
      check for a reaction file name of "." or "./",
      if so these will flag as an error.
    */
    if (output_filename_length < 3) {
      if ((strcmp(output_filename,"./") == 0)  || (strcmp(output_filename,".") == 0) || (output_filename_length < 1)) {
	success = 0;
	fprintf(stderr,"create_output_filenames: Error invalid rxn_file was %s\n",output_filename);
	fflush(stderr);
      }
    }
  }
  if (success) {
    if ((output_filename_base_length + 11) > state->max_filename_len) {
      output_filename_base_length = state->max_filename_len - 12;
      /*
	Probably should print a warning mewssage here.
      */
      fprintf(stderr,"create_output_filenames: Warning, truncated output_filename_base to %ld characters\n",output_filename_base_length);
      fflush(stderr);
    }
    state->output_filename_base_length = output_filename_base_length;
    if (init_conc_filename[0] == '\0') {
      strncpy(init_conc_filename,output_filename,output_filename_base_length);
      strcpy((char*)&init_conc_filename[output_filename_base_length],".concs");
    }
    if (log_filename[0] == '\0') {
      strncpy(log_filename,output_filename,output_filename_base_length);
      strcpy((char*)&log_filename[output_filename_base_length],".log");
    }
    if (output_filename_empty) {
      /*
	No output file was specified so we use the rxn_file prefix as the
	output file base, and make its suffix .out.
      */
      strcpy((char*)&output_filename[output_filename_base_length],".out");
    }
    if (counts_out_filename[0] == '\0') {
      strncpy(counts_out_filename,output_filename,output_filename_base_length);
      strcpy((char*)&counts_out_filename[output_filename_base_length],".count");
    }
    if (concs_out_filename[0] == '\0') {
      strncpy(concs_out_filename,output_filename,output_filename_base_length);
      strcpy((char*)&concs_out_filename[output_filename_base_length],".concs");
    }
    if (ode_concs_filename[0] == '\0') {
      strncpy(ode_concs_filename,output_filename,output_filename_base_length);      
      strcpy((char*)&ode_concs_filename[output_filename_base_length],".ode_concs");
    }
    if (ode_counts_filename[0] == '\0') {
      strncpy(ode_counts_filename,output_filename,output_filename_base_length);      
      strcpy((char*)&ode_counts_filename[output_filename_base_length],".ode_counts");
    }
    strncpy(ode_dconcs_filename,output_filename,output_filename_base_length);
    strcpy((char*)&ode_dconcs_filename[output_filename_base_length],".ode_dconcs");
    strncpy(ode_lklhd_filename,output_filename,output_filename_base_length);
    strcpy((char*)&ode_lklhd_filename[output_filename_base_length],".ode_lklhd");
    strncpy(ode_bflux_filename,output_filename,output_filename_base_length);
    strcpy((char*)&ode_bflux_filename[output_filename_base_length],".ode_bflux");
    strncpy(net_lklhd_filename,output_filename,output_filename_base_length);
    strcpy((char*)&net_lklhd_filename[output_filename_base_length],".nlklhd");
    strncpy(nl_bndry_flx_filename,output_filename,output_filename_base_length);
    strcpy((char*)&nl_bndry_flx_filename[output_filename_base_length],".nl_flux");
    strncpy(aux_data_filename,output_filename,output_filename_base_length);
    strcpy((char*)&aux_data_filename[output_filename_base_length],".aux");
    if (rxn_lklhd_filename[0] == '\0') {
      strncpy(rxn_lklhd_filename,output_filename,output_filename_base_length);
      strcpy((char*)&rxn_lklhd_filename[output_filename_base_length],".lklhd");
    }
    if (free_energy_filename[0] == '\0') {
      strncpy(free_energy_filename,output_filename,output_filename_base_length);
      strcpy((char*)&free_energy_filename[output_filename_base_length],".fe");
    }
    if (restart_filename[0] == '\0') {
      strncpy(restart_filename,output_filename,output_filename_base_length);
      strcpy((char*)&restart_filename[output_filename_base_length],".rstrt");
    }
    if (rxn_view_filename[0] == '\0') {
      strncpy(rxn_view_filename,output_filename,output_filename_base_length);
      strcpy((char*)&rxn_view_filename[output_filename_base_length],".view");
    }
    if (bndry_flux_filename[0] == '\0') {
      strncpy(bndry_flux_filename,output_filename,output_filename_base_length);
      strcpy((char*)&bndry_flux_filename[output_filename_base_length],".flux");
    }
    if (rxn_echo_filename[0] == '\0') {
      strncpy(rxn_echo_filename,output_filename,output_filename_base_length);
      strcpy((char*)&rxn_echo_filename[output_filename_base_length],".echo");
    }
    if (rxn_mat_filename[0] == '\0') {
      strncpy(rxn_mat_filename,output_filename,output_filename_base_length);
      strcpy((char*)&rxn_mat_filename[output_filename_base_length],".mat");
    }
    if (dg0ke_filename[0] == '\0') {
      strncpy(dg0ke_filename,output_filename,output_filename_base_length);
      strcpy((char*)&dg0ke_filename[output_filename_base_length],".dg0ke");
    }
    if (dictionary_filename[0] == '\0') {
      strncpy(dictionary_filename,output_filename,output_filename_base_length);
      strcpy((char*)&dictionary_filename[output_filename_base_length],".dict");
    }
    if (ode_sens_filename[0] == '\0') {
      strncpy(ode_sens_filename,output_filename,output_filename_base_length);
      strcpy((char*)&ode_sens_filename[output_filename_base_length],".sens");
    }
    if (ode_dsens_filename[0] == '\0') {
      strncpy(ode_dsens_filename,output_filename,output_filename_base_length);
      strcpy((char*)&ode_dsens_filename[output_filename_base_length],".dsens");
    }
  }
  return(success);
}
