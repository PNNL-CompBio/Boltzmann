/* parse_reactions_file.c
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
#include "count_ws.h"
#include "count_nws.h"
#include "upcase.h"
#include "is_a_coef.h"

#include "parse_reactions_file.h"
int parse_reactions_file(struct state_struct *state) {
  /*
    Determine the number of reactions,
    Total length of molecules names, and
    max possiblie number of molecules, and
    total length of the file, total length
    of compartment names, total length of pathway_names,
    total_length reaction titles.
    Called by: boltzmann_init
    Calls    : parse_rxn_file_keyword
               count_ws
	       count_nws
	       upcase
	       fseek, feof, fgets, fprintf, fflush, sscanf (intrinsic)
  */
  struct rxn_struct *reactions;
  struct rxn_struct *reaction;
  struct rxn_matrix_struct *rxns_matrix;
  struct istring_elem_struct *unsorted_molecules;
  struct istring_elem_struct *rxn_molecules;
  struct istring_elem_struct *unsorted_cmpts;
  int64_t *keyword_lens;
  int64_t *rxn_ptrs;
  int64_t *molecules_indices;
  int64_t *coefficients;
  int64_t rxn_buff_len;
  int64_t total_length;
  int64_t molecules_len;
  int64_t rxn_title_len;
  int64_t pathway_len;
  int64_t compartment_len;
  int64_t line_len;
  int64_t seek_offset;
  int64_t rxn_title_pos;
  int64_t pathway_pos;
  int64_t compartment_pos;
  int64_t molecules_pos;
  int64_t align_len; 
  int64_t align_mask; 
  int64_t len;
  int64_t sll;
  char **matrix_text;
  char *rxn_buffer;
  char *fgp;
  char **keywords;
  char *rctnts;
  char *prdcts;
  char *rxn_title_text;
  char *pathway_text;
  char *compartment_text;
  char *molecules_text;
  char *raw_molecules_text;

  int success;
  int rxns;

  int molecules;
  int ws_chars;

  int line_type;
  int kl;

  int cmpts;
  int pos;

  int skip;
  int sl;

  int ns;
  int j;

  int mol_pos;
  int mol_pos_lim;

  int padding;
  int ci;

  int colon_loc;
  int padi;

  FILE *rxn_fp;
  FILE *lfp;
  success = 1;
  seek_offset = (int64_t)0;
  rxn_buff_len = state->rxn_buff_len;
  lfp          = state->lfp;
  rxn_fp       = state->rxn_fp;
  align_len    = state->align_len;
  align_mask   = state->align_mask;
  if (rxn_fp == NULL) {
    success = 0;
    fprintf(stderr,"parse_reactions_file: reaction file not open, %s\n",
	    state->reaction_file);
    fflush(stderr);
  }
  if (success) {
    /*
      Seek to beginning of file.
    */
    fseek(rxn_fp,seek_offset,SEEK_SET);
    rxn_buffer = state->rxn_buffer;
    keywords   = state->rxn_file_keywords;
    keyword_lens = state->rxn_file_keyword_lengths;
    rxns    = 0;
    molecules = 0;
    cmpts   = 0;
    /*
      Should get a reaction line, a pathway line, a left line, a right line,
      a dgzero line, a dgzero-units line and a terminating // line per
      reaction. We build a state machine
      This is not quite accurate,as some reactions may not
      have a pathway line, and Bill wants to also add an additional
      compartment line.
      
    */
    rxn_title_pos          = (int64_t)0;
    pathway_pos            = (int64_t)0;
    compartment_pos        = (int64_t)0;
    molecules_pos            = (int64_t)0;
    unsorted_molecules            = state->unsorted_molecules;
    unsorted_cmpts              = state->unsorted_cmpts;
    rxn_title_text              = state->rxn_title_text;
    pathway_text                = state->pathway_text;
    compartment_text            = state->compartment_text;
    molecules_text                = state->molecules_text;
    raw_molecules_text            = state->raw_molecules_text;
    reactions                   = state->reactions;
    rxns_matrix                 = state->reactions_matrix;
    rxn_ptrs                    = rxns_matrix->rxn_ptrs;
    molecules_indices              = rxns_matrix->molecules_indices;
    coefficients                = rxns_matrix->coefficients;
    matrix_text                 = rxns_matrix->text;
    reaction                    = reactions;
    reaction->lcompartment      = NULL;
    reaction->rcompartment      = NULL;
    reaction->pathway           = NULL;
    reaction->left_compartment  = -1;
    reaction->right_compartment = -1;
    reaction->num_reactants     = 0;
    reaction->num_products      = 0;
    rxn_ptrs[rxns]              = molecules;
    fgp = fgets(rxn_buffer,rxn_buff_len,rxn_fp);
    state->max_molecule_len = 0;
    state->min_molecule_len = rxn_buff_len;
    state->max_compartment_len = 0;
    state->min_compartment_len = rxn_buff_len;
    while ((fgp && success) && (! feof(rxn_fp))) {
      line_len = strlen(rxn_buffer);
      /*
	Check that last character in line is a newline and replace with 
	a \0.
      */
      if (rxn_buffer[line_len-1] == '\n') {
	rxn_buffer[line_len-1] = '\0';
	line_len -= 1;
      } else {
	fprintf(stderr,"parse_reactions_file: Error input line longer than"
		" %d characters\n",rxn_buff_len);
	fflush(stderr);
	success = 0;
	break;
      }
      if (success) {
      	/*
      	  state, every reaction starts with a REACTION line
      	  ends with a // line.
      	  It requires: LEFT, RIGHT, DGZERO and DGZERO-UNITS lines.
      	  It may have PATHWAY and LEFT_COMPARTMENT, RIGHT_COMPARTMENT  lines.

      	  Allow pre REACTION lines in reactions file as a header.
      	  So before the first REACTION line the state is -1
      		
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
	    Get size of reaction title allow +1 for terminating \0.
          */
	  rxn_title_len = (int64_t)(line_len - kl - ws_chars+1);
	  /*
	    Copy the reaction title.
	  */
	  reaction->title = (char *)&rxn_title_text[rxn_title_pos];
	  strcpy(reaction->title,(char*)&rxn_buffer[kl+ws_chars]);
	  /*
	    caution bit twiddle follows: 
	    padding = (align_len - (rxn_titl_len % align_len)) % align_len
	    Thou shalt not use %.
	  */
	  padding = (align_len - (rxn_title_len & align_mask)) & align_mask;
	  rxn_title_pos += rxn_title_len + padding;
	  break;
	case 1:
	  /*
	    Pathway line.
	  */
	  pathway_len = line_len - kl - ws_chars + 1;
	  reaction->pathway = (char *)&pathway_text[pathway_pos];
	  strcpy(reaction->pathway,(char*)&rxn_buffer[kl+ws_chars]);
	  padding = (align_len - (pathway_len & align_mask)) & align_mask;
	  pathway_pos += pathway_len + padding;
	  break;
	case 2: 
	  /*
	    Compartment line.
	  */
	  compartment_len = line_len - kl - ws_chars + 1;
	  if (compartment_len > state->max_compartment_len) {
	    state->max_compartment_len = compartment_len;
	  } else {
	    if (compartment_len < state->min_compartment_len) {
	      state->min_compartment_len = compartment_len;
	    }
	  }
	  reaction->lcompartment = (char *)&compartment_text[compartment_pos];
	  strcpy(reaction->lcompartment,(char*)&rxn_buffer[kl+ws_chars]);
	  upcase(compartment_len,reaction->lcompartment,
		 reaction->lcompartment);
	  if ((strcmp(reaction->lcompartment,"V") == 0) ||
	      (strcmp(reaction->lcompartment,"C") == 0)) {
	    fprintf(stderr,"parse_reactions_file: Error compartments may not have single character names V or C\n");
	    fflush(stderr);
	    success = 0;
	    break;
	  }
	  reaction->rcompartment = NULL;
	  padding = (align_len - (compartment_len & align_mask)) & align_mask;
	  compartment_pos += compartment_len + padding;
	  reaction->left_compartment = cmpts;
	  reaction->right_compartment = cmpts;
	  unsorted_cmpts->string = reaction->lcompartment;
	  unsorted_cmpts->m_index  = -1;
	  unsorted_cmpts->c_index  = cmpts;
	  unsorted_cmpts += 1; /* Caution address arithmetic */
	  cmpts += 1;
	  break;
	case 3: 
	  /*
	    Left Compartment line.
	  */
	  compartment_len = line_len - kl - ws_chars + 1;
	  if (compartment_len > state->max_compartment_len) {
	    state->max_compartment_len = compartment_len;
	  } else {
	    if (compartment_len < state->min_compartment_len) {
	      state->min_compartment_len = compartment_len;
	    }
	  }
	  reaction->lcompartment = (char *)&compartment_text[compartment_pos];
	  strcpy(reaction->lcompartment,(char*)&rxn_buffer[kl+ws_chars]);
	  upcase(compartment_len,reaction->lcompartment,
		 reaction->lcompartment);
	  padding = (align_len - (compartment_len & align_mask)) & align_mask;
	  compartment_pos += compartment_len + padding;
	  reaction->left_compartment = cmpts;
	  unsorted_cmpts->string = reaction->lcompartment;
	  unsorted_cmpts->c_index  = cmpts;
	  unsorted_cmpts += 1; /* Caution address arithmetic */
	  cmpts += 1;
	  break;
	case 4: 
	  /*
	    Right Compartment line.
	  */
	  compartment_len = line_len - kl - ws_chars + 1;
	  if (compartment_len > state->max_compartment_len) {
	    state->max_compartment_len = compartment_len;
	  } else {
	    if (compartment_len < state->min_compartment_len) {
	      state->min_compartment_len = compartment_len;
	    }
	  }
	  reaction->rcompartment = (char *)&compartment_text[compartment_pos];
	  strcpy(reaction->rcompartment,(char*)&rxn_buffer[kl+ws_chars]);
	  upcase(compartment_len,reaction->rcompartment,
		 reaction->rcompartment);
	  padding = (align_len - (compartment_len & align_mask)) & align_mask;
	  compartment_pos += compartment_len + padding;
	  reaction->right_compartment = cmpts;
	  unsorted_cmpts->string = reaction->rcompartment;
	  unsorted_cmpts->c_index  = cmpts;
	  unsorted_cmpts += 1; /* Caution address arithmetic */
	  cmpts += 1;
	  break;
	case 5:
	  /*
	    A left line, count and record reactant molecules and coefficients.
	  */
	  rctnts = (char *)&rxn_buffer[kl+ ws_chars];
	  pos = 0;
	  len = strlen(rctnts);
	  while (pos < len) {
	    sl = (int64_t)count_nws(rctnts);
	    if (sl > 0) {
	      molecules_indices[molecules] = molecules;
	      if (is_a_coef(sl,rctnts)) {
		rctnts[sl] = '\0';
		coefficients[molecules] = - atoi(rctnts);
		rctnts[sl] = ' ';

		pos += sl;
		rctnts += sl; /* Caution address arithmetic. */

		skip   =  count_ws(rctnts);

		pos += skip;
		rctnts += skip; /* Caution address arithmetic. */

		sl = count_nws(rctnts);
	      } else {
		coefficients[molecules] = -1;
	      }
	      /*
		At this juncture we need to check for
		a "local" compartment specification of the form 
		:compartment trailing the molecule name. We do not
		allow spaces on either side of the semicolon.
	      */
	      ci = -1;
	      rctnts[sl] = '\0';
	      colon_loc = find_colon(rctnts);
	      if (colon_loc >= 0) {
		/* 
		  We had a local :compartment attached.
		  determine its length, store it and shorten the
		  length for the molecule - do not forget to 
		  allow for the terminating null.
		*/
		compartment_len = sl - colon_loc;
		sl = colon_loc;
		rctnts[sl] = '\0';
		if (compartment_len > state->max_compartment_len) {
		  state->max_compartment_len = compartment_len;
		} else {
		  if (compartment_len < state->min_compartment_len) {
		    state->min_compartment_len = compartment_len;
		  }
		}
		unsorted_cmpts->string = (char *)&compartment_text[compartment_pos];
		unsorted_cmpts->c_index  = cmpts;
		unsorted_cmpts += 1; /* Caution address arithmetic */
		ci = cmpts;
		cmpts += 1;
		strcpy(unsorted_cmpts->string,(char*)&rctnts[colon_loc+1]);
		upcase(compartment_len,unsorted_cmpts->string,
		       unsorted_cmpts->string);
		padding = (align_len - (compartment_len & align_mask)) & 
		  align_mask;
		compartment_pos += compartment_len + padding;
	      }
	      if (sl > state->max_molecule_len) {
		state->max_molecule_len = sl;
	      } else {
		if (sl < state->min_molecule_len) {
		  state->min_molecule_len = sl;
		}
	      }
	      matrix_text[molecules] = (char*)&raw_molecules_text[molecules_pos];
	      sll = (int64_t)sl + (int64_t)1;
	      padding = (align_len - (sll & align_mask)) & align_mask;
	      strcpy((char *)&raw_molecules_text[molecules_pos],rctnts);
	      upcase(sl,(char *)&raw_molecules_text[molecules_pos],
		     (char *)&molecules_text[molecules_pos]);
	      unsorted_molecules->string = (char *)&molecules_text[molecules_pos];
	      unsorted_molecules->m_index  = molecules;
	      unsorted_molecules->c_index  = ci;
	      unsorted_molecules += 1; /* Caution address arithmetic. */

	      molecules_pos += (int64_t)(sll + padding);
	      rctnts[sl] = ' ';

	      pos += sl; /* Caution address arithmetic. */
	      rctnts += sl;
	      molecules += 1;
	      reaction->num_reactants += 1;
	      skip = count_ws(rctnts);

	      pos += skip; /* Caution address arithmetic. */
	      rctnts += skip;

	      if (pos < len) {
		if (rctnts[0] != '+') {
		  fprintf(stderr,"parse_reactions_file: Error, missing + "
			  "between reactants, reactions file line was \n%s\n",
			  rxn_buffer);
		  fflush(stderr);
		  success = 0;
		  break;
		}

		pos += 1;
		rctnts += 1; /* Caution address arithmetic.*/

		skip = count_ws(rctnts);

		pos  += skip;
		rctnts += skip; /* Caution address arithmetic.*/
	      }
	    } /* end if (sl > 0)  - a reactant was found . */
	  } /* end while (pos > len) */
	  break;
	case 6:
	  /*
	    A right line, count product molecules.
	  */
	  prdcts = (char *)&rxn_buffer[kl + ws_chars];
	  pos = 0;
	  len = strlen(prdcts);
	  while (pos < len) {
	    sl = (int64_t)count_nws(prdcts);
	    if (sl > 0) {
	      molecules_indices[molecules] = molecules;
	      if (is_a_coef(sl,prdcts)) {
		prdcts[sl] = '\0';
		coefficients[molecules] = atoi(prdcts);
		prdcts[sl] = ' ';

		pos += sl;
		prdcts += sl; /* Caution address arithmetic. */

		skip   =  count_ws(prdcts);
		
		pos += skip;
		prdcts += skip; /* Caution address arithmetic.*/

		sl = count_nws(prdcts);
	      } else {
		coefficients[molecules] = 1;
	      }
	      ci = -1;
	      prdcts[sl] = '\0';
	      colon_loc = find_colon(prdcts);
	      if (colon_loc >= 0) {
		/* 
		  We had a local :compartment attached.
		  determine its length, store it and shorten the
		  length for the molecule - do not forget to 
		  allow for the terminating null.
		*/
		compartment_len = sl - colon_loc;
		sl = colon_loc;
		rctnts[sl] = '\0';
		if (compartment_len > state->max_compartment_len) {
		  state->max_compartment_len = compartment_len;
		} else {
		  if (compartment_len < state->min_compartment_len) {
		    state->min_compartment_len = compartment_len;
		  }
		}
		unsorted_cmpts->string = (char *)&compartment_text[compartment_pos];
		unsorted_cmpts->c_index  = cmpts;
		unsorted_cmpts += 1; /* Caution address arithmetic */
		ci = cmpts;
		cmpts += 1;
		strcpy(unsorted_cmpts->string,(char*)&rctnts[colon_loc+1]);
		upcase(compartment_len,unsorted_cmpts->string,
		       unsorted_cmpts->string);
		padding = (align_len - (compartment_len & align_mask)) & 
		  align_mask;
		compartment_pos += compartment_len + padding;
	      }
	      if (sl > state->max_molecule_len) {
		state->max_molecule_len = sl;
	      } else {
		if (sl < state->min_molecule_len) {
		  state->min_molecule_len = sl;
		}
	      }
	      if (sl > state->max_molecule_len) {
		state->max_molecule_len = sl;
	      } else {
		if (sl < state->min_molecule_len) {
		  state->min_molecule_len = sl;
		}
	      }
	      matrix_text[molecules] = (char*)&raw_molecules_text[molecules_pos];
	      prdcts[sl] = '\0';
	      sll = (int64_t)sl + (int64_t)1;
	      padding = (align_len - (sll & align_mask)) & align_mask;
	      strcpy((char *)&raw_molecules_text[molecules_pos],prdcts);
	      upcase(sl,(char *)&raw_molecules_text[molecules_pos],
		     (char *)&molecules_text[molecules_pos]);
	      unsorted_molecules->string = (char *)&molecules_text[molecules_pos];
	      unsorted_molecules->m_index  = molecules;
	      unsorted_molecules->c_index  = ci;
	      unsorted_molecules += 1; /* Caution address arithmetic. */
	      molecules_pos += (int64_t)(sll + padding);
	      prdcts[sl] = ' ';

	      pos += sl;
	      prdcts += sl; /* Caution address arithmetic. */

	      molecules += 1;
	      reaction->num_products += 1;
	      skip = count_ws(prdcts);

	      pos += skip;
	      prdcts += skip; /* Caution address arithmetic. */

	      if (pos < len) {
		if (prdcts[0] != '+') {
		  fprintf(stderr,"parse_reactions_file: Error, missing + "
			  "between products, reactions file line was \n%s\n",
			  rxn_buffer);
		  fflush(stderr);
		  success = 0;
		  break;
		}
		pos    += 1;
		prdcts += 1; /* Caution address arithmetic. */

		skip = count_ws(prdcts);

		pos  += skip;
		prdcts += skip; /* Caution address arithmetic. */
	      }
	    } /* end if (sl > 0)  - a reactant was found . */
	  } /* end while (pos > len) */
	  break;
	case 7:
	  /* 
	     A DGZERO line
	  */
	  ns = sscanf ((char*)&rxn_buffer[ws_chars+kl],"%le",
		       &reaction->delta_g0);
	  if (ns < 1) {
	    fprintf(stderr,
		    "parse_reactions_file: malformed DGZERO line was\n%s\n",
		    rxn_buffer);
	    fflush(stderr);
	    success = 0;
	    break;
	  }
	  break;
	case 8:
	  /* 
	     A DGZERO-UNITS line
	  */
	  sl = count_nws((char*)&rxn_buffer[ws_chars+kl]);
	  upcase(sl,(char*)&rxn_buffer[ws_chars+kl],(char*)&rxn_buffer[ws_chars+kl]);
	  if (strncmp((char*)&rxn_buffer[ws_chars+kl],"KJ/MOL",6) == 0) {
	    reaction->unit_i = 1;
	  } else {
	    reaction->unit_i = 0;
	  }
	  break;
	case 9:
	  reaction->self_id = rxns;
	  /*
	    Since a compartment line could have come any where in
	    the reaction input lines, we need to go back
	    and properly set the compartments for each of the
	    molecules in the reaction. Reactants get the
	    reaction->left_compartment value (for c_index) and
	    products get the reaction->right_compartment value for c_index;
	  */
	  mol_pos = rxn_ptrs[rxns];
	  rxn_molecules = (struct istring_elem_struct *)&state->unsorted_molecules[mol_pos];
	  mol_pos_lim = mol_pos + reaction->num_reactants +
	    reaction->num_products;
	  for (j = mol_pos;j<mol_pos_lim;j++) {
	    /*
	      Only look to set the compartment index if the molecule did
	      not have a local compartment (:compartment).
	    */
	    if (rxn_molecules->c_index < 0) {
	      if (coefficients[j] < 0) {
		rxn_molecules->c_index = reaction->left_compartment;
	      } else {
		rxn_molecules->c_index = reaction->right_compartment;
	      }
	    }
	    rxn_molecules += 1; /* Caution address arithmetic */
	  }
	  rxns += 1;
	  /*
	    Caution address arithmetic follows
	    reaction = (struct rxn_struct*)&reactions[rxns];
	  */
	  reaction += 1;
	  if (rxns < state->number_reactions) {
	    reaction->lcompartment      = NULL;
	    reaction->rcompartment      = NULL;
	    reaction->pathway           = NULL;
	    reaction->left_compartment  = -1;
	    reaction->right_compartment = -1;
	    reaction->num_reactants     = 0;
	    reaction->num_products      = 0;
	    rxn_ptrs[rxns]              = molecules;
	  }
	  break;
        default:
	break;
	}
	fgp = fgets(rxn_buffer,rxn_buff_len,rxn_fp);
      }/* end if (success) */
    } /* end while(fgp...) */
  } /* end if (success) */
  rxn_ptrs[rxns] = molecules;
  /*
    Check that last line was a //.
  */
  if (line_type != state->num_rxn_file_keywords - 1) {
    fprintf(stderr,
	    "parse_reactions_file: Error reactions file did not end in //\n");
    fflush(stderr);
    success = 0;
  }
  /*
    These are now set in size_rxns_file.
  state->reaction_file_length = total_length;
  state->number_reactions = rxns;
  state->number_molecules   = molecules;
  state->number_compartments = cmpts;
  state->molecules_len = molecules_len;
  state->pathway_len = pathway_len;
  state->compartment_len = compartment_len;
  state->rxn_title_len = rxn_title_len;
  */
  return(success);
}
