/* init_rxn_file_keywords.c
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

#include "init_rxn_file_keywords.h"
void init_rxn_file_keywords(struct state_struct *state) {
  /*
    Initialize the reaction file keywords and lengths.
    Called by: size_rxns_file.
    Calls:     strncpy(intrinsic)
  */
  char *keyword_buffer;
  char **keywords;
  int64_t *keyword_len;

  keyword_buffer = state->rxn_file_keyword_buffer;
  keywords       = state->rxn_file_keywords;
  keyword_len    = state->rxn_file_keyword_lengths;
  strcpy(keyword_buffer,
	 "REACTIONPATHWAY COMPARTMENT     LEFT_COMPARTMENT");
  strcpy((char *)&keyword_buffer[48],
	 "RIGHT_COMPARTMENT       LEFT    RIGHT   DGZERO  ");
  strcpy((char *)&keyword_buffer[96],
	 "DGZERO-UNITS    ACTIVITY        PREGULATION     ");
  strcpy((char *)&keyword_buffer[144],
	 "NREGULATION     ENZYME_LEVEL    //              ");
  keywords[0]  	     = (char*)(&keyword_buffer[0]);
  keywords[1]  	     = (char*)(&keyword_buffer[8]);
  keywords[2]  	     = (char*)(&keyword_buffer[16]);
  keywords[3]  	     = (char*)(&keyword_buffer[32]);
  keywords[4]  	     = (char*)(&keyword_buffer[48]);
  keywords[5]  	     = (char*)(&keyword_buffer[72]);
  keywords[6]  	     = (char*)(&keyword_buffer[80]);
  keywords[7]  	     = (char*)(&keyword_buffer[88]);
  keywords[8]  	     = (char*)(&keyword_buffer[96]);
  keywords[9]  	     = (char*)(&keyword_buffer[112]);
  keywords[10] 	     = (char*)(&keyword_buffer[128]);
  keywords[11] 	     = (char*)(&keyword_buffer[144]);
  keywords[12] 	     = (char*)(&keyword_buffer[160]);
  keywords[13] 	     = (char*)(&keyword_buffer[176]);
  keyword_len[0]     = (int64_t)8;
  keyword_len[1]     = (int64_t)7;
  keyword_len[2]     = (int64_t)11;
  keyword_len[3]     = (int64_t)16;
  keyword_len[4]     = (int64_t)17;
  keyword_len[5]     = (int64_t)4;
  keyword_len[6]     = (int64_t)5;
  keyword_len[7]     = (int64_t)6;
  keyword_len[8]     = (int64_t)12;
  keyword_len[9]     = (int64_t)8;
  keyword_len[10]    = (int64_t)11;
  keyword_len[11]    = (int64_t)11;
  keyword_len[12]    = (int64_t)12;
  keyword_len[13]    = (int64_t)2;
}
