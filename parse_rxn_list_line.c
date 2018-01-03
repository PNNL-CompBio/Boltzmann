#include "boltzmann_structs.h"

#include "count_ws.h"
#include "count_nws.h"

#include "parse_rxn_list_line.h"
int parse_rxn_list_line(struct state_struct *local_state, struct boot_state_struct *boot_state, int64_t rxn_file) {
  /*
    Parse a line from a reaction list file to set the 
    reaction_file, init_conc_file,  and compartment_file or sbml_file fields
    of local_state. Workspace is provided by boot_state.
    
    Called by: boltzmann_boot
    Calls:     count_ws,count_nws,fgets,strncpy
    
  */
  int64_t rxn_list_line_len;
  int success;
  int skip1;

  int skip2;
  int skip3;

  int rxn_file_name_len;
  int concs_file_name_len;
  int cmpts_file_name_len;
  int padi;

  char *rp;
  char *rxn_list_buff;
  char *reaction_file;
  char *init_conc_file;
  char *compartment_file;
  char *sbml_file;
  char *rxn_file_name;
  char *concs_file_name;
  char *cmpts_file_name;
  
  FILE *rxn_list_fp;
  FILE *lfp;

  success = 1;
  lfp               = boot_state->lfp;
  rxn_list_fp       = boot_state->rxn_list_fp;
  rxn_list_buff     = boot_state->rxn_list_buffer;
  rxn_list_line_len = boot_state->rxn_list_buffer_len;
  reaction_file     = local_state->reaction_file;
  init_conc_file    = local_state->init_conc_file;
  compartment_file  = local_state->compartment_file;
  sbml_file         = local_state->sbml_file;
  rp = fgets(rxn_list_buff,rxn_list_line_len,rxn_list_fp);
  if (rp == NULL) {
    success = 0;
    if (lfp) {
      fprintf(lfp,"parse_rxn_list_line: Error reading %lld'th line of rxns_list file\n",
	      rxn_file);
      fflush(lfp);
    }
  } else {
    skip1 = count_ws(rxn_list_buff);
    rxn_file_name = &rxn_list_buff[skip1];
    rxn_file_name_len = count_nws(rxn_file_name);
    strncpy(reaction_file,rxn_file_name,rxn_file_name_len);
    reaction_file[rxn_file_name_len] = '\0';
    skip2 = count_ws(&rxn_file_name[rxn_file_name_len]);
    concs_file_name = &rxn_file_name[rxn_file_name_len+skip2];
    concs_file_name_len = count_nws(concs_file_name);
    if (concs_file_name_len > 0) {
      sbml_file[0] = '\0';
      strncpy(init_conc_file,concs_file_name,concs_file_name_len);
      init_conc_file[concs_file_name_len] = '\0';
      skip3 = count_ws(&concs_file_name[concs_file_name_len]);
      cmpts_file_name = &concs_file_name[concs_file_name_len+ skip3];
      cmpts_file_name_len = count_nws(cmpts_file_name);
      if (cmpts_file_name_len > 0) {
	strncpy(compartment_file,cmpts_file_name,cmpts_file_name_len);
	compartment_file[cmpts_file_name_len] = '\0';
      }
    } else {
      /*
	Only one file name was given we will assume it was an
	SBML file 
      */
      strncpy(sbml_file,rxn_file_name,rxn_file_name_len);
      sbml_file[rxn_file_name_len] = '\0';
    }
  }
  return(success);
}

