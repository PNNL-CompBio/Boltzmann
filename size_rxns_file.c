/* size_rxns_file.c
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
#include "count_species.h"
#include "count_ws.h"

#include "size_rxns_file.h"
int size_rxns_file(struct state_struct *state) {
  /*
    Determine the number of reactions,
    Total length of species names, and
    max possiblie number of species, and
    total length of the file. Total length
    of compartment names.
    Called by: boltzmann
    Calls    : count_ws
  */
  int64_t rxn_buff_len;
  int64_t total_length;
  int64_t species_len;
  int64_t rxn_title_len;
  int64_t pathway_len;
  int64_t compartment_len;
  int64_t line_len;
  int64_t *keyword_lens;
  char *rxn_buffer;
  char *fgp;
  char **keywords;
  char *header[8];
  char *rctnts;
  char *prdcts;

  int success;
  int rxns;

  int species;
  int ws_chars;

  int line_type;
  int kl;

  int cmpts;
  int pad1;

  FILE *rxn_fp;
  FILE *lfp;
  success = 1;
  rxn_buff_len = state->rxn_buff_len;
  lfp          = state->lfp;
  rxn_fp       = fopen(state->reaction_file,"r");
  if (rxn_fp == NULL) {
    success = 0;
    fprintf(stderr,"size_rxn_file: unable to open reaction file, %s\n",
	    state->reaction_file);
    fflush(stderr);
  }
  if (success) {
    rxn_buffer = state->rxn_buffer;
    init_rxn_file_keywords(state);
    keywords   = state->rxn_file_keywords;
    keyword_lens = state->rxn_file_keyword_lengths;
    rxns = 0;
    species = 0;
    cmpts   = 0;
    /*
      Should get a reaction line, a pathway line, a left line, a right line,
      a dgzero line, a dgzero-units line and a terminating // line per
      reaction. We build a state machine
      This is not quite accurate,as some reactions may not
      have a pathway line, and Bill wants to also add an additional
      compartment line.
      
    */
    total_length = (int64_t)0;
    species_len  = (int64_t)0;
    pathway_len  = (int64_t)0;
    compartment_len = (int64_t)0;
    rxn_title_len  = (int64_t)0;
    fgp = fgets(rxn_buffer,rxn_buff_len,rxn_fp);
    while (fgp && (! feof(rxn_fp))) {
      line_len     =  strlen(rxn_buffer);
      total_length += line_len;
      /*
	state, every reaction starts with a REACTION line
	ends with a // line.
	It requires: LEFT, RIGHT, DGZERO and DGZERO-UNITS lines.
	It may have PATHWAY and COMPARTMENT lines.

	Allow pre REACTION lines in reactions file as a header.
	So before the first REACTION line the state is -1

	When a reaction line is seen then the state should be 0
	and any of the other lines may be seen. Need to record
	which types have been seen, when parsing the rxns file.
	For sizing as we are doing here, it doesn't matter.
	
	When a // is seen then we hit transfer state 1, where 
	we must get a REACTION line next.
      */
      line_type = parse_rxn_file_keyword(rxn_buffer,state);
      if (line_type >= 0) {
	kl = keyword_lens[line_type];
	ws_chars = count_ws((char *)&rxn_buffer[kl]);      

      } else {
	kl       = 0;
	ws_chars = 0;
      }
      /*
	Line header matched.
      */
      switch (line_type) {
	case 0: 
	  /*
	    Reaction title line.
	    Get size of reaction title.
          */
	  rxns += 1;
	  rxn_title_len += line_len - kl - ws_chars;
	  break;
	case 1:
	  /*
	    Pathway line.
	  */
	  ws_chars = count_ws((char *)&rxn_buffer[kl]);
	  pathway_len += line_len - kl - ws_chars;
	  break;
	case 2: 
	  /*
	    Compartment line.
	  */
	  cmpts += 1;
	  compartment_len += line_len - kl - ws_chars;
	  break;
	case 3:
	  /*
	    A left line, count species.
	  */
	  rctnts = (char *)&rxn_buffer[kl];
	  species += count_species(rctnts,&species_len);
	  break;
	case 4:
	  /*
	    A right line, count species.
	  */
	  prdcts = (char *)&rxn_buffer[kl];
	  species += count_species(prdcts,&species_len);
	  break;
	case 5:
	case 6:
	case 7:
        default:
	  break;
      }
      fgp = fgets(rxn_buffer,rxn_buff_len,rxn_fp);
    } /* end while(fgp...) */
  } /* end if (success) */
  state->reaction_file_length = total_length;
  state->number_reactions = rxns;
  state->number_species   = species;
  state->number_compartments = cmpts;
  state->species_len = species_len;
  state->pathway_len = pathway_len;
  state->compartment_len = compartment_len;
  state->rxn_title_len = rxn_title_len;
  return(success);
}
