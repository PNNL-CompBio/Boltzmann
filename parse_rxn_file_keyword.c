/* parse_rxn_file_keyword.c
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

#include "parse_rxn_file_keyword.h"

int parse_rxn_file_keyword(char *rxn_buffer,struct state_struct *state){
  /*
    Parse the reaction file keyword theoretically at the start of 
    rxn_buffer.
    Called by: size_rxns_file
  */
  char **keywords;
  char *keyword;
  int64_t *keyword_len;
  int i;
  int num_rxn_file_keywords;
  int line_type;
  int kl;
  int skip;
  int ki;
  FILE *lfp;
  line_type = -1;
  keywords = state->rxn_file_keywords;
  keyword_len = state->rxn_file_keyword_lengths;
  num_rxn_file_keywords = state->num_rxn_file_keywords;
  skip = count_ws(rxn_buffer);
  keyword = (char *)&rxn_buffer[skip];
  kl   = count_nws(keyword);
  for (i=0;i<num_rxn_file_keywords;i++) {
    ki = keyword_len[i];
    if (ki == kl) {
      if (strncmp(keyword,keywords[i],ki) == 0) {
	line_type = i;
	break;
      }
    }
  }
  if (line_type < 0) {
    lfp = state->lfp;
    if (lfp) {
      fprintf(lfp,"parse_rxn_file_keyword: Unrecognized keyword\n%s\n",
	      rxn_buffer);
    }
  }
  return(line_type);
}
