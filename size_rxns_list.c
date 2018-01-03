/* size_rxns_list.c
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

#include "size_rxns_list.h"
int size_rxns_list(struct boot_state_struct *boot_state) {
  /*
    Determine the number of reaction files,
    Called by: boot_init
    Calls    : fopen, fgets, fprintf, fflush 
  */
  int64_t rxn_buff_len;
  char    *rxn_list_buffer;
  char    *rxn_list_file;
  char    *fgp;

  int success;
  int rxns;

  int num_reaction_files;
  int padi;

  FILE *rxn_fp;
  FILE *lfp;

  success = 1;
  lfp             = boot_state->lfp;
  rxn_buff_len    = boot_state->rxn_list_buffer_len;
  rxn_list_buffer = boot_state->rxn_list_buffer;
  rxn_list_file   = boot_state->rxn_list_file;
  rxn_fp          = fopen(rxn_list_file,"r");
  if (rxn_fp == NULL) {
    success = 0;
    num_reaction_files = -1;
    if (lfp) {
      fprintf(lfp,"size_rxn_list: unable to open reaction list file, %s\n",
	      rxn_list_file);
      fflush(lfp);
    }
  }
  if (success) {
    num_reaction_files = 0;
    while (! feof(rxn_fp)) {
      fgp = fgets(rxn_list_buffer,rxn_buff_len,rxn_fp);
      if (fgp) {
	num_reaction_files += 1;
      }
    }
  }
  fclose(rxn_fp);
  boot_state->num_reaction_files = num_reaction_files;
  return(success);
}
