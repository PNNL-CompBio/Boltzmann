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

#include "boltzmann_structs.h"

#include "init_rxn_file_keywords.h"
#include "parse_rxn_file_keyword.h"
#include "count_molecules_and_cmpts.h"
#include "count_ws.h"
#include "count_nws.h"

#include "size_rxns_file.h"
int size_rxns_file(struct state_struct *state,
		   char *reaction_file) {
  /*
    Determine the number of reactions,
    Total length of molecules names, and
    max possible number of molecules, and
    total length of the file. Total length
    of compartment names.
    Called by: io_size_init
    Calls    : init_rxn_file_keywords
               parse_rxn_file_keyword
               count_ws,
	       count_nws,
               count_molecules_and_cmpts,
               fopen, fgets, fprintf, fflush (intrinsic)
  */
  int64_t rxn_buff_len;
  int64_t total_length;
  int64_t molecules_len;
  int64_t rxn_title_len;
  int64_t pathway_len;
  int64_t compartment_len;
  int64_t line_len;
  int64_t align_len;
  int64_t max_regs_per_rxn;
  int64_t regulation_len;
  int64_t *keyword_lens;
  char *rxn_buffer;
  char *fgp;
  char **keywords;
  char *rctnts;
  char *prdcts;

  int success;
  int rxns;

  int molecules;
  int species_len;

  int ws_chars;
  int line_type;

  int kl;
  int cmpts;

  int line_no;
  int padi;

  FILE *rxn_fp;
  FILE *lfp;
  success = 1;
  rxn_buff_len = state->max_param_line_len << 1;
  align_len    = state->align_len;
  lfp          = state->lfp;
  rxn_fp       = fopen(reaction_file,"r");
  if (rxn_fp == NULL) {
    success = 0;
    fprintf(stderr,"size_rxn_file: unable to open reaction file, %s\n",
	    reaction_file);
    fflush(stderr);
  }
  if (success) {
    /*
    state->rxn_fp = rxn_fp;
    */
    rxn_buffer = state->param_buffer;
    init_rxn_file_keywords(state);
    keywords   = state->rxn_file_keywords;
    keyword_lens = state->rxn_file_keyword_lengths;
    rxns = 0;
    molecules = 0;
    cmpts   = 0;
    /*
      Should get a reaction line, a pathway line, a left line, a right line,
      a dgzero line, a dgzero-units line and a terminating // line per
      reaction. We build a state machine
      This is not quite accurate,as some reactions may not
      have a pathway line, and Bill wants to also add an additional
      compartment, left_compartment, right_compartement  lines.
      
    */
    total_length = (int64_t)0;
    molecules_len  = (int64_t)0;
    pathway_len  = (int64_t)0;
    regulation_len = (int64_t)0;
    /*
      Allow for space for the empty compartment null string. 
    */
    compartment_len = (int64_t)align_len;
    if (compartment_len < 1) {
      compartment_len = 1;
    }
    rxn_title_len  = (int64_t)0;
    fgp = fgets(rxn_buffer,rxn_buff_len,rxn_fp);
    line_no = 1;
    while (fgp && (! feof(rxn_fp))) {
      line_len     =  strlen(rxn_buffer);
      total_length += line_len;
      /*
	state, every reaction starts with a REACTION line
	ends with a // line.
	It requires: LEFT, RIGHT, DGZERO and DGZERO-UNITS lines.
	It may have PATHWAY and COMPARTMENT lines,
	and PREGULATION and NREGULATION and ACTIVITY lines.

	Allow pre REACTION lines in reactions file as a header.
	So before the first REACTION line the state is -1

	When a reaction line is seen then the state should be 0
	and any of the other lines may be seen. Need to record
	which types have been seen, when parsing the rxns file.
	For sizing as we are doing here, it doesn't matter.
	
	When a // is seen then we hit transfer state 1, where 
	we must get a REACTION line next.
      */
      line_type = parse_rxn_file_keyword(rxn_buffer,line_no,state);
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
      case 0: /* REACTION */
	  /*
	    Reaction title line.
	    Get size of reaction title.
          */
	  rxns += 1;
	  rxn_title_len += line_len - kl - ws_chars;
	  break;
      case 1: /* PATHWAY */
	  /*
	    Pathway line.
	  */
	  pathway_len += line_len - kl - ws_chars;
	  break;
      case 2: /* COMPARTMENT */ 
      case 3: /* LEFT_COMPARTMENT */
      case 4: /* RIGHT_COMPARTMENT */
	  /*
	    Compartment line.
	  */
	  cmpts += 1;
	  compartment_len += line_len - kl - ws_chars;
	  break;
      case 5: /* LEFT */
	  /*
	    A left line, count molecules.
	  */
   	  rctnts = (char *)&rxn_buffer[kl];
	  count_molecules_and_cmpts(rctnts,&molecules,&cmpts,&molecules_len,
				    &compartment_len);
	  /*
	    Here we also need to count the number of compartments
	    (number of : on the line.
	  */
	  break;
      case 6: /* RIGHT */
	  /*
	    A right line, count molecules.
	  */
	  prdcts = (char *)&rxn_buffer[kl];
	  count_molecules_and_cmpts(prdcts,&molecules,&cmpts,&molecules_len,
				    &compartment_len);
	  /*
	    Here we also need to count the number of compartments
	    (number of : on the line.
	  */
	  break;
      case 7: /* DGZERO */
      case 8: /* DGZERO-UNITS */
	  break;
      case 9:  /* ACTIVITY */
      case 12: /* ENZYME_LEVEL */
      case 13: /* // */
      case 14: /* FORWARD_RATE */
      case 15: /* REVERSE_RATE */
	  break;
      case 10: /* PREGULATION */
      case 11: /* NREGLATION */
	  /*
	    A PREGULATION or NREGULATION line.
	  */
	  species_len = count_nws((char*)&rxn_buffer[kl+ws_chars]);
	  regulation_len += species_len;
      default:
	break;
      }
      fgp = fgets(rxn_buffer,rxn_buff_len,rxn_fp);
      line_no += 1;
    } /* end while(fgp...) */
    state->reaction_file_length = total_length;
    state->number_reactions = rxns;
    state->number_molecules   = molecules;
    /*
      Always have an empty compartment.
    */
    cmpts += 1;
    /*
      We need to add space for padding done in parse_rxns
    */
    align_len = state->align_len;
    max_regs_per_rxn = state->max_regs_per_rxn;
    molecules_len += molecules * align_len;
    compartment_len += cmpts * align_len;
    pathway_len += rxns*align_len;
    rxn_title_len += rxns * align_len;
    regulation_len += rxns * max_regs_per_rxn * align_len;
    state->number_compartments = cmpts;
    state->molecule_text_length    = molecules_len;
    state->pathway_text_length     = pathway_len;
    state->compartment_text_length = compartment_len;
    state->reaction_titles_length   = rxn_title_len;
    state->regulation_text_length  = regulation_len;
    fclose(rxn_fp);
  }
  return(success);
}
