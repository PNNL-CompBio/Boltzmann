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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <float.h>
#include <signal.h>
#include <unistd.h>

#include "boltzmann_structs.h"

#include "init_rxn_file_keywords.h"
#include "parse_rxn_file_keyword.h"
#include "count_molecules.h"
#include "count_ws.h"

#include "size_rxns_list.h"
int size_rxns_list(struct state_struct *state) {
  /*
    Determine the number of reaction files,
    Called by: boltzmann_boot
    Calls    : count_ws,
               count_molecules,
               fopen, fgets, fprintf, fflush (intrinsic)
  */
  int64_t rxn_buff_len;
  char *rxn_buffer;
  char *fgp;

  int success;
  int rxns;

  int num_reaction_files;
  int padi;

  FILE *rxn_fp;
  FILE *lfp;
  success = 1;
  rxn_buff_len = state->max_param_line_len << 1;
  lfp          = state->lfp;
  rxn_fp       = fopen(state->reaction_file,"r");
  if (rxn_fp == NULL) {
    success = 0;
    num_reaction_files = -1;
    fprintf(stderr,"size_rxn_list: unable to open reaction list file, %s\n",
	    state->reaction_file);
    fflush(stderr);
  }
  if (success) {
    rxn_buffer = state->param_buffer;
    num_reaction_files = 0;
    while (! feof(rxn_fp)) {
      fgp = fgets(rxn_buffer,rxn_buff_len,rxn_fp);
      if (fgp) {
	num_reaction_files += 1;
      }
    }
  }
  fclose(rxn_fp);
  return(num_reaction_files);
}
