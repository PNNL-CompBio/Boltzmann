/* sbml_set_file_names.c
*******************************************************************************
boltzmann

Pacific Northwest National Laboratory, Richland, WA 99352.

Copyright (c) 2010 Battelle Memorial Institute.

Publications based on work performed using the software should include 
the following citation as a reference:


Licensed under the Educational Community License, Version 2.0 (the "License"); 
you may not use this file except in compliance with the License. 
The terms and conditions of the License may be found in 
ECL-2.0_LICENSE_TERMS.TXT in the directory containing this file.
        
Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR 
CONDITIONS OF ANY KIND, either express or implied. See the License for the 
specific language governing permissions and limitations under the License.
******************************************************************************/
#include "boltzmann_structs.h"
#include "sbml_set_file_names.h"

int sbml_set_file_names(struct sbml2bo_struct *state) {
  /*
    Set the file names for the rxns.dat, concs.in and cmpts.dat files.
    Assumes sbml file name is in sbml_file field of state.
    Called by sbml_to_boltzmann and sbml2bo
  */
  int success;
  int base_len;
  int in_len;
  int max_file_name_len;
  char *sbml_file;
  char *concs_in_file;
  char *rxns_dat_file;
  char *cmpts_dat_file;
  char *log_file;
  char *tail;
  sbml_file      = state->sbml_file;
  concs_in_file  = state->concs_in_file;
  rxns_dat_file  = state->rxns_dat_file;
  cmpts_dat_file = state->cmpts_dat_file;
  log_file       = state->log_file;
  max_file_name_len = state->file_name_len;
  /* Determine the length of the base name - this is just the sbml file name
     minus the .sbml if it ends in one, else the full file name.
  */
  if (success) {
    in_len = strlen(sbml_file);
    tail = (char*)&sbml_file[in_len-5];
    if (strncmp(tail,".sbml",5) == 0) {
      base_len = in_len - 5;
    } else {
      base_len = in_len;
    }
    /*
      Use an _concs.in suffix for the concentrations input file
      an _rxns.dat suffix for the reactions data file,
      an _cmpts.dat suffix for the compartments data file.
      and a .log suffix for the log file.
    */
    if (base_len + 10 > max_file_name_len) {
      base_len = max_file_name_len - 10;
    }
    if (base_len < 0) {
      success = 0;
    }
  }
  if (success) {
    strncpy(concs_in_file,sbml_file,base_len);
    tail = (char*)&concs_in_file[base_len];
    strcpy (tail,"_concs.in");
    strncpy(rxns_dat_file,sbml_file,base_len);
    tail = (char*)&rxns_dat_file[base_len];
    strcpy (tail,"_rxns.dat");
    strncpy(cmpts_dat_file,sbml_file,base_len);
    tail = (char*)&cmpts_dat_file[base_len];
    strcpy (tail,"_cmpts.dat");
    strncpy(log_file,sbml_file,base_len);
    tail = (char*)&log_file[base_len];
    strcpy (tail,".log");
  }
  return(success);
}
