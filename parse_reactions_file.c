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
#include "boltzmann_structs.h"

#include "parse_rxn_file_keyword.h"
#include "count_ws.h"
#include "count_nws.h"
#include "upcase.h"
#include "is_a_coef.h"
#include "parse_side_line.h"

#include "parse_reactions_file.h"
int parse_reactions_file(struct state_struct *state,
			 char *reaction_file) {
  /*

    This routine fills the unsorted_molecules and unsorted_cmpts
    molecule structures. 
    It fills the reactions structure.
    It also sets the molecules_indices, coefficients,
    and text fields of the reactions matrix.
    It also fills the rxn_title_text, pathway_text, compartment_text
    and molecules_text buffers.
   
    An molecule_struct has 4 fields. The first three are set by 
    this routine.
    
        string  - pointer to a null terminated string that
	          contains the molecule or compartment string.
	m_index - If string points to a molecule, m_index
         	  is the ordinal position of the molecule in the
		  LEFT and RIGHT lines of the reactions.dat input file.
		  If the string is a compartment name m_index will be -1.
        c_index - If the string points to a compartment, c_index is the
	          ordinal position of the compartment name in the
		  COMPARTMENT, LEFT_COMPARTMENT, RIGHT_COMPARTMENT line, or
		  colon preceded fields of the LEFT and RIGHT lines in 
		  the reactions.dat input file.
		  If the string points to a molecule and the molecule
		  is in a compartment, c_index is set to the ordinal number
		  of the containing compartment, otherwise if the 
		  molecule is not in a compartment its c_index is -1.

	variable - An indicator set by the read_initial_concentrations
	          routine indicating whether or not this molecule is 
		  held fixed in concentration.

    Called by: rxns_init
    Calls    : parse_rxn_file_keyword,
               count_ws,
	       count_nws,
	       is_a_coef,
	       upcase,
	       fseek, feof, fgets, fprintf, fflush, sscanf (intrinsic)
  */
  struct reaction_struct *reactions;
  struct reaction_struct *reaction;
  struct reactions_matrix_struct *rxns_matrix;
  struct molecule_struct *unsorted_molecules;
  struct molecule_struct *rxn_molecules;
  struct compartment_struct *unsorted_cmpts;
  double *activities;
  double *reg_constant;
  double *reg_exponent;
  double *reg_drctn;
  int64_t *keyword_lens;
  int64_t *rxn_ptrs;
  int64_t *molecules_indices;
  int64_t *coefficients;
  int64_t *reg_species;
  int64_t rxn_buff_len;
  int64_t rxn_title_len;
  int64_t pathway_len;
  int64_t compartment_len;
  int64_t line_len;
  int64_t seek_offset;
  int64_t rxn_title_pos;
  int64_t pathway_pos;
  int64_t compartment_pos;
  int64_t molecules_pos;
  int64_t regulation_pos;
  int64_t align_len; 
  int64_t align_mask; 
  int64_t reg_base;
  int64_t nr;
  char *rxn_buffer;
  char *fgp;
  char **keywords;
  char *rctnts;
  char *prdcts;
  char *rxn_title_text;
  char *pathway_text;
  char *compartment_text;
  char *regulation_text;
  char *molecules_text;
  char *raw_molecules_text;
  char *title;
  char *lcompartment;
  char *rcompartment;
  char *pathway;
  char *metabolite;

  size_t word1_len;

  int success;
  int rxns;

  int molecules;
  int regulations;

  int ws_chars;
  int line_type;

  int kl;
  int cmpts;

  int sl;
  int ns;

  int j;
  int mol_pos;

  int mol_pos_lim;
  int padding;

  int side;
  int word1;

  int i;
  int max_regs_per_rxn;

  int reg_count;
  int reg_pos;


  FILE *rxn_fp;
  FILE *lfp;
  success = 1;
  seek_offset  	   = (int64_t)0;
  rxn_buff_len 	   = state->max_param_line_len<<1;
  lfp          	   = state->lfp;
  align_len    	   = state->align_len;
  align_mask   	   = state->align_mask;
  reg_constant 	   = state->reg_constant;
  reg_exponent 	   = state->reg_exponent;
  reg_species  	   = state->reg_species;
  reg_drctn    	   = state->reg_drctn;
  max_regs_per_rxn = (int)state->max_regs_per_rxn;
  reg_base         = (int64_t)0;
  reg_pos          = reg_base;
  /*
  strcpy(state->reaction_file,reaction_file);
  */
  rxn_fp           = fopen(reaction_file,"r");
  if (rxn_fp == NULL) {
    success = 0;
    if (lfp) {
      fprintf(lfp,"parse_reactions_file: reaction file not open, %s\n",
	      state->reaction_file);
      fflush(lfp);
    }
  }
  if (success) {
    /*
      Seek to beginning of file.
    fseek(rxn_fp,seek_offset,SEEK_SET);
    */
    rxn_buffer = state->param_buffer;
    keywords   = state->rxn_file_keywords;
    keyword_lens = state->rxn_file_keyword_lengths;
    rxns    = 0;
    molecules = 0;
    regulations = 0;
    /*
      Should get a reaction line, a pathway line, a left line, a right line,
      a dgzero line, a dgzero-units line and a terminating // line per
      reaction. We build a state machine
      This is not quite accurate,as some reactions may not
      have a pathway line, and Bill wants to also add an additional
      compartment line.
      
    */
    rxn_title_pos               = (int64_t)0;
    pathway_pos                 = (int64_t)0;
    compartment_pos             = (int64_t)0;
    molecules_pos               = (int64_t)0;
    regulation_pos              = (int64_t)0;
    
    unsorted_molecules          = state->unsorted_molecules;
    unsorted_cmpts              = state->unsorted_cmpts;
    rxn_title_text              = state->rxn_title_text;
    pathway_text                = state->pathway_text;
    compartment_text            = state->compartment_text;
    molecules_text              = state->molecules_text;
    raw_molecules_text          = state->raw_molecules_text;
    regulation_text             = state->regulation_text;
    activities                  = state->activities;
    reactions                   = state->reactions;
    rxns_matrix                 = state->reactions_matrix;
    rxn_ptrs                    = rxns_matrix->rxn_ptrs;
    molecules_indices           = rxns_matrix->molecules_indices;
    coefficients                = rxns_matrix->coefficients;
    reaction                    = reactions;
    reaction->lcompartment      = 0;
    reaction->rcompartment      = 0;
    reaction->pathway           = -1;
    reaction->left_compartment  = 0;
    reaction->right_compartment = 0;
    reaction->num_reactants     = 0;
    reaction->num_products      = 0;
    reaction->activity          = 1.0;
    reaction->temp_kelvin       = state->temp_kelvin;
    reaction->deltag0_computed  = 0;
    reaction->ph                = state->ph;
    reaction->ionic_strength    = state->ionic_strength;
    activities[0]               = 1.0;
    for (i=0;i<max_regs_per_rxn;i++) {
      reg_constant[reg_base+i] = 1.0;
      reg_exponent[reg_base+i] = 0.0;
      reg_drctn[reg_base+i]    = 1.0;
      reg_species[reg_base+i]  = (int64_t)(-1);
    }
    reg_count                  = 0;
    rxn_ptrs[rxns]              = molecules;
    fgp = fgets(rxn_buffer,rxn_buff_len,rxn_fp);
    state->max_molecule_len = (int64_t)0;
    state->min_molecule_len = rxn_buff_len;
    state->max_compartment_len = (int64_t)0;
    state->min_compartment_len = rxn_buff_len;
    /*
      Build in the empty compartment.
    */
    unsorted_cmpts->string   = compartment_pos;
    unsorted_cmpts->c_index  = 0;
    unsorted_cmpts           += 1; /* Caution address arithmetic */
    compartment_text[0]      = '\0';
    regulation_text[0]       = '\0';
    compartment_pos          = align_len;
    cmpts = 1;
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
	success = 0;
	if (lfp) {
	  fprintf(lfp,"parse_reactions_file: Error input line longer than"
		" %lld characters\n",rxn_buff_len);
	  fflush(lfp);
	}
	break;
      }
      if (success) {
      	/*
      	  state, every reaction starts with a REACTION line
      	  ends with a // line.
      	  It requires: LEFT, RIGHT, DGZERO and DGZERO-UNITS lines.
      	  It may have PATHWAY and LEFT_COMPARTMENT, RIGHT_COMPARTMENT  
	  ACTIVITY, PREGULATION and or NREGULATION lines.

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
	word1 = kl + ws_chars;
	word1_len = (size_t)count_nws((char *)&rxn_buffer[word1]);
	/*
	  Line header matched.
	*/
	switch (line_type) {
	case 0: 
	  /*
	    Reaction title line.
	    Get size of reaction title allow +1 for terminating \0.
          */
	  rxn_title_len = (int64_t)(line_len - word1 + 1);
	  /*
	    Copy the reaction title.
	  reaction->title = (char *)&rxn_title_text[rxn_title_pos];
	  */
	  reaction->title = rxn_title_pos;
	  title  = (char *)&rxn_title_text[rxn_title_pos];
	  strncpy(title,(char*)&rxn_buffer[word1],word1_len);
	  title[word1_len] = '\0';
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
	  pathway_len = line_len - word1 + 1;
	  /*
	  reaction->pathway = (char *)&pathway_text[pathway_pos];
	  */
	  reaction->pathway = pathway_pos;
	  pathway = (char *)&pathway_text[pathway_pos];
	  strncpy(pathway,(char*)&rxn_buffer[word1],word1_len);
	  pathway[word1_len] = '\0';
	  padding = (align_len - (pathway_len & align_mask)) & align_mask;
	  pathway_pos += pathway_len + padding;
	  break;
	case 2: 
	  /*
	    Compartment line.
	  */
	  compartment_len = word1_len;
	  if ((int64_t)compartment_len > state->max_compartment_len) {
	    state->max_compartment_len = (int64_t)compartment_len;
	  } else {
	    if ((int64_t)compartment_len < state->min_compartment_len) {
	      state->min_compartment_len = (int64_t)compartment_len;
	    }
	  }
	  /*
	  reaction->lcompartment = (char *)&compartment_text[compartment_pos];
	  */
	  reaction->lcompartment = compartment_pos;
	  lcompartment = (char *)&compartment_text[compartment_pos];
	  strncpy(lcompartment,(char*)&rxn_buffer[word1],word1_len);
	  lcompartment[word1_len] = '\0';
	  upcase(compartment_len,lcompartment,
		 lcompartment);
	  if ((strcmp(lcompartment,"V") == 0) ||
	      (strcmp(lcompartment,"C") == 0)) {
	    if (lfp) {
	      fprintf(lfp,"parse_reactions_file: Error compartments may not have single character names V or C\n");
	      fflush(lfp);
	    }
	    success = 0;
	    break;
	  }
	  reaction->rcompartment = 0;
	  padding = (align_len - (compartment_len & align_mask)) & align_mask;
	  reaction->left_compartment = cmpts;
	  reaction->right_compartment = cmpts;
	  /*
	  unsorted_cmpts->string = lcompartment;
	  */
	  unsorted_cmpts->string = compartment_pos;
	  unsorted_cmpts->volume = 0.0;
	  unsorted_cmpts->c_index  = cmpts;
	  unsorted_cmpts += 1; /* Caution address arithmetic */
	  compartment_pos += compartment_len + padding;
	  cmpts += 1;
	  break;
	case 3: 
	  /*
	    Left Compartment line.
	  */
	  compartment_len = word1_len;
	  if ((int64_t)compartment_len > state->max_compartment_len) {
	    state->max_compartment_len = (int64_t)compartment_len;
	  } else {
	    if ((int64_t)compartment_len < state->min_compartment_len) {
	      state->min_compartment_len = (int64_t)compartment_len;
	    }
	  }
	  /*
	  reaction->lcompartment = (char *)&compartment_text[compartment_pos];
	  */
	  reaction->lcompartment = compartment_pos;
	  lcompartment = (char *)&compartment_text[compartment_pos];
	  strncpy(lcompartment,(char*)&rxn_buffer[word1],word1_len);
	  lcompartment[word1_len] = '\0';
	  upcase(compartment_len,lcompartment,
		 lcompartment);
	  padding = (align_len - (compartment_len & align_mask)) & align_mask;
	  reaction->left_compartment = cmpts;
	  /*
	  unsorted_cmpts->string = lcompartment;
	  */
	  unsorted_cmpts->string = compartment_pos;
	  unsorted_cmpts->volume = 0.0;
	  unsorted_cmpts->c_index  = cmpts;
	  unsorted_cmpts += 1; /* Caution address arithmetic */
	  compartment_pos += compartment_len + padding;
	  cmpts += 1;
	  break;
	case 4: 
	  /*
	    Right Compartment line.
	  */
	  compartment_len = word1_len;
	  if ((int64_t)compartment_len > state->max_compartment_len) {
	    state->max_compartment_len = (int64_t)compartment_len;
	  } else {
	    if ((int64_t)compartment_len < state->min_compartment_len) {
	      state->min_compartment_len = (int64_t)compartment_len;
	    }
	  }
	  /*
	  reaction->rcompartment = (char *)&compartment_text[compartment_pos];
	  */
	  reaction->rcompartment = compartment_pos;
	  rcompartment = (char *)&compartment_text[compartment_pos];
	  strncpy(rcompartment,(char*)&rxn_buffer[word1],word1_len);
	  rcompartment[word1_len] = '\0';
	  upcase(compartment_len,rcompartment,
		 rcompartment);
	  padding = (align_len - (compartment_len & align_mask)) & align_mask;
	  reaction->right_compartment = cmpts;
	  /*
	  unsorted_cmpts->string = rcompartment;
	  */
	  unsorted_cmpts->string = compartment_pos;
	  unsorted_cmpts->volume = 0.0;
	  unsorted_cmpts->c_index  = cmpts;
	  unsorted_cmpts += 1; /* Caution address arithmetic */
	  compartment_pos += compartment_len + padding;
	  cmpts += 1;
	  break;
	case 5:
	  /*
	    A left line, count and record reactant molecules and coefficients.
	  */
	  side = -1;
	  rctnts = (char *)&rxn_buffer[word1];
	  success = parse_side_line(rctnts,(int64_t *)&molecules_pos,
				    (int64_t *)&compartment_pos,
				    (int *)&molecules,
				    (int *)&cmpts,
				    (struct reaction_struct *)reaction,
				    state,
				    side);
	  break;
	case 6:
	  /*
	    A right line, count product molecules.
	  */
	  prdcts = (char *)&rxn_buffer[word1];
	  side   = 1;
	  success = parse_side_line(prdcts,(int64_t *)&molecules_pos,
				    (int64_t *)&compartment_pos,
				    (int *)&molecules,
				    (int *)&cmpts,
				    (struct reaction_struct *)reaction,
				    state,
				    side);
	  break;
	case 7:
	  /* 
	     A DGZERO line
	  */
	  ns = sscanf ((char*)&rxn_buffer[word1],"%le",
		       &reaction->delta_g0);
	  if (ns < 1) {
	    title  = (char *)&rxn_title_text[reaction->title];
	    if (lfp) {
	      fprintf(lfp,
		      "parse_reactions_file: Error: malformed DGZERO line"
		      " for reaction %s was\n%s\n",
		      title,rxn_buffer);
	      fflush(lfp);
	    }
	    success = 0;
	    break;
	  }
	  break;
	case 8:
	  /* 
	     A DGZERO-UNITS line
	  */
	  sl = word1_len;
	  if (sl < 1) {
	    title  = (char *)&rxn_title_text[reaction->title];
	    if (lfp) {
	      fprintf(lfp,
		    "parse_reactions_file: Error: malformed DGZERO-UNITS line,"
		    " for reaction %s, was %s, using KJ/MOL\n",
		    title,rxn_buffer);
	      fflush(lfp);
	    }
	    reaction->unit_i = 1;
	  } else {
	    upcase(sl,(char*)&rxn_buffer[word1],(char*)&rxn_buffer[word1]);
	    if (strncmp((char*)&rxn_buffer[word1],"KJ/MOL",6) == 0) {
	      reaction->unit_i = 1;
	    } else {
	      reaction->unit_i = 0;
	    }
	  }
	  break;
        case 9:
	  /*
	    An ACTIVITY line.
	  */
	  ns = sscanf ((char*)&rxn_buffer[word1],"%le",
		       &reaction->activity);
	  if (ns < 1) {
	    if (lfp) {
	      fprintf(lfp,
		    "parse_reactions_file: malformed ACTIVITY line was\n%s\n",
		    rxn_buffer);
	      fflush(lfp);
	    }
	    success = 0;
	    break;
	  }
	  break;
	case 10:
	  /*
	    PREGULATION LINE
	  */
	  if (reg_count >= max_regs_per_rxn - 1) {
	    success = 0;
	    if (lfp) {
	      fprintf(lfp,"parse_reactions_file: Error limit of %d regulations per reaction was exceeded. Increase the MAX_REGS_PER_RXN parameter and rerun.\n",
		      max_regs_per_rxn);
	      fflush(lfp);
	    }
	    break;
	  }
	  metabolite = (char *)&regulation_text[regulation_pos];
	  strncpy(metabolite,(char*)&rxn_buffer[word1],word1_len);
	  metabolite[word1_len] = '\0';
	  upcase(word1_len,metabolite,metabolite);
	  padding = (align_len - (word1_len & align_mask)) & align_mask;
	  /*
	    reg_pos = reg_base + reg_count;
	  */
	  reg_species[reg_pos] = regulation_pos;
	  reg_drctn[reg_pos]   = 1.0;
	  nr = sscanf((char*)&rxn_buffer[word1+word1_len],"%le %le",&reg_constant[reg_pos],&reg_exponent[reg_pos]);
	  if (nr < 2) {
	    success = 0;
	    if (lfp) {
	      fprintf(lfp,"parse_reactions_file: Error PREGULATION LINE, %s\n"
		      "did not have valid constant and exponent.\n",
		      rxn_buffer);
	      fflush(lfp);
	    }
	    break;
	  }
	  regulation_pos += (word1_len + padding);
	  reg_count += 1;
	  reg_pos += 1;		 
	  break;
	case 11:
	  /*
	    NREGULATION LINE
	  */
	  if (reg_count >= max_regs_per_rxn - 1) {
	    success = 0;
	    if (lfp) {
	      fprintf(lfp,"parse_reactions_file: Error limit of %d regulations per reaction was exceeded. Increase the MAX_REGS_PER_RXN parameter and rerun.\n",
		      max_regs_per_rxn);
	      fflush(lfp);
	    }
	    break;
	  }
	  metabolite = (char *)&regulation_text[regulation_pos];
	  strncpy(metabolite,(char*)&rxn_buffer[word1],word1_len);
	  metabolite[word1_len] = '\0';
	  upcase(word1_len,metabolite,metabolite);
	  padding = (align_len - (word1_len & align_mask)) & align_mask;
	  /*
	    reg_pos = reg_base + reg_count;
	  */
	  reg_species[reg_pos] = regulation_pos;
	  reg_drctn[reg_pos]   = 0.0;
	  nr = sscanf((char*)&rxn_buffer[word1+word1_len],"%le %le",&reg_constant[reg_pos],&reg_exponent[reg_pos]);
	  if (nr < 2) {
	    success = 0;
	    if (lfp) {
	      fprintf(lfp,"parse_reactions_file: Error NREGULATION LINE, %s\n"
		      "did not have valid constant and exponent.\n",
		      rxn_buffer);
	      fflush(lfp);
	    }
	    break;
	  }
	  regulation_pos += (word1_len + padding);
	  reg_count += 1;
	  reg_pos += 1;		 

	  break;
	case 12:
	  /*
	    // reaction terminator line
	  */
	  reaction->self_id = rxns;
	  activities[rxns]  = reaction->activity;
	  /*
	    Since a compartment line could have come any where in
	    the reaction input lines, we need to go back
	    and properly set the compartments for each of the
	    molecules in the reaction. Reactants get the
	    reaction->left_compartment value (for c_index) and
	    products get the reaction->right_compartment value for c_index;
	  */
	  mol_pos = rxn_ptrs[rxns];
	  rxn_molecules = (struct molecule_struct *)&state->unsorted_molecules[mol_pos];
	  mol_pos_lim = mol_pos + reaction->num_reactants +
	    reaction->num_products;
	  for (j = mol_pos;j<mol_pos_lim;j++) {
	    /*
	      Only look to set the compartment index if the molecule did
	      not have a local compartment (:compartment).
	    */
	    if (rxn_molecules->c_index == 0) {
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
	    reaction = (struct reaction_struct*)&reactions[rxns];
	  */
	  reaction += 1;
	  reg_base  += max_regs_per_rxn;
	  reg_pos   =  reg_base;
	  reg_count = 0;
	  if (rxns < (int)state->number_reactions) {
	    reaction->lcompartment      = 0;
	    reaction->rcompartment      = 0;
	    reaction->pathway           = -1;
	    reaction->left_compartment  = 0;
	    reaction->right_compartment = 0;
	    reaction->num_reactants     = 0;
	    reaction->num_products      = 0;
	    reaction->activity          = 1.0;
	    /*
	      The following three lines added by DGT on 4/18/2013
	    */
	    reaction->temp_kelvin       = state->temp_kelvin;
	    reaction->ph                = state->ph;
	    reaction->ionic_strength    = state->ionic_strength;
	    reaction->deltag0_computed  = 0;

	    rxn_ptrs[rxns]              = molecules;
	    for (i=0;i<max_regs_per_rxn;i++) {
	      reg_constant[reg_base+i] = 1.0;
	      reg_exponent[reg_base+i] = 0.0;
	      reg_drctn[reg_base+i]    = 1.0;
	      reg_species[reg_base+i]  = (int64_t)(-1);
	    }
	  }
	  break;
        default:
	break;
	} /* end switch(line_type) */
      }/* end if (success) */
      fgp = fgets(rxn_buffer,rxn_buff_len,rxn_fp);
    } /* end while(fgp...) */
    rxn_ptrs[rxns] = molecules;
    /*
      Check that last line was a //.
    */
    if (success) {
      if (line_type != (int)(state->num_rxn_file_keywords) - 1) {
	if (lfp) {
	  fprintf(lfp,
		  "parse_reactions_file: Error reactions file did not end in //\n");
	  fflush(lfp);
	}
	success = 0;
      }
    } else {
      if (lfp) {
	fprintf(lfp,"parse_reactions_file: Error occured reading reactions file %s, see above.\n",reaction_file);
	fflush(lfp);
      }
    }
    if (rxn_fp != NULL) {
      fclose(rxn_fp);
    }
  } /* end if (success) */
  return(success);
}
