/* parse_reactions_file.c
*******************************************************************************
boltzmann

Pacific Northwest National Laboratory, Richland, WA 99352.

Copyright (c) 2018 Battelle Memorial Institute.

******************************************************************************/
#include "boltzmann_structs.h"

#include "parse_rxn_file_keyword.h"
#include "count_ws.h"
#include "count_nws.h"
#include "upcase.h"
#include "parse_side_line.h"
#include "boltzmann_compress_reactions.h"

#include "parse_reactions_file.h"
int parse_reactions_file(struct state_struct *state,
			 char *reaction_file) {
  /*

    This routine fills the unsorted_molecules and unsorted_cmpts
    molecule structures. For the compartments it sets the volumes
    to be 0 - to indicate they have not been set.
    Then the volume, recip_volume, ph and ionic_strength fields of
    the compartment are set in species init in read_compartment_sizes,
    and read_initial_concentrations.
    It fills the reactions structure, and builds the reactions matrix.
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
               parse_side_line
	       upcase,
	       fseek, feof, fgets, fprintf, fflush, sscanf, strncpy,
	       strcmp, strncmp
	       
    Uses the following fields in state
      lfp,
      max_param_line_len,
      max_regs_per_rxn
      align_mask,
      align_len,

    Sets the following fields in state:
      reg_constant,
      reg_exponent,
      reg_species,
      reg_drctn,
      coeff_sum,
      use_rxn,
      reactions,
      reactions_matrix, and subfields
      pathway_text,
      compartment_text,
      molecules_text,
      raw_molecules_text,
      regulation_text
      activities,
      enzyme_level
     
   and the unsorted_molecules struct, the rxn_molecules struct, and 
   the usnorted_compartments fiedls 
  */
  struct reaction_struct *reactions;
  struct reaction_struct *reaction;
  struct reactions_matrix_struct *rxns_matrix;
  struct molecule_struct *unsorted_molecules;
  struct molecule_struct *rxn_molecules;
  struct compartment_struct *unsorted_cmpts;
  struct compartment_struct *compartment;
  double *activities;
  double *enzyme_level;
  double *reg_constant;
  double *reg_exponent;
  double *reg_drctn;
  double *forward_rc;
  double *reverse_rc;
  double *recip_coeffs;
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
  int64_t max_rxn_title_pos;
  int64_t max_pathway_pos;
  int64_t max_compartment_pos;
  int64_t max_regulation_pos;
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
  int  *coeff_sum;
  int  *use_rxn;

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

  int last_line_type;
  int sum_coeffs;

  int rxn_use;
  int rxn_not_used;

  int num_products;
  int num_reactants;

  int had_a_title;
  int had_a_dg0;

  int had_units;
  int line_no;

  int max_cmpt;
  int rtpd;

  int ppd;
  int cpd;

  int rgpd;
  int padi;
  
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
  coeff_sum        = state->coeff_sum;
  use_rxn          = state->use_rxn;
  max_cmpt         = state->number_compartments - 1;
  max_rxn_title_pos = state->reaction_titles_length;
  max_pathway_pos   = state->pathway_text_length;
  max_compartment_pos = state->compartment_text_length;
  max_regulation_pos  = state->regulation_text_length;
  max_regs_per_rxn = (int)state->max_regs_per_rxn;
  reg_base         = (int64_t)0;
  reg_pos          = reg_base;
  line_type        = 0;
  rxn_not_used     = 0;

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
      We add optional K_FORWARD and K_REVERSE lines for use 
      with DELTA_CONCS_CHOICE 9 in the ode parts.
      
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
    enzyme_level                = state->enzyme_level;
    reactions                   = state->reactions;
    rxns_matrix                 = state->reactions_matrix;
    forward_rc                  = state->forward_rc;
    reverse_rc                  = state->reverse_rc;
    rxn_ptrs                    = rxns_matrix->rxn_ptrs;
    molecules_indices           = rxns_matrix->molecules_indices;
    coefficients                = rxns_matrix->coefficients;
    recip_coeffs                = rxns_matrix->recip_coeffs;
    reaction                    = reactions;
    /*
      Initialize first reaction.
    */
    reaction->lcompartment      = 0;
    reaction->rcompartment      = 0;
    reaction->pathway           = -1;
    reaction->left_compartment  = 0;
    reaction->right_compartment = 0;
    reaction->num_reactants     = 0;
    reaction->num_products      = 0;
    reaction->activity          = 1.0;
    reaction->enzyme_level      = 1.0;
    reaction->temp_kelvin       = state->temp_kelvin;
    reaction->ph                = state->ph;
    reaction->ionic_strength    = state->ionic_strength;
    reaction->forward_rc        = -1.0;
    reaction->reverse_rc        = -1.0;
    reaction->unit_i            = 1;
    reaction->delta_g0          = 0.0;
    /*
      We will change semantics here. reaction->deltag0_computed 
      will be initialized to -1. 
      If a reaction has a DGZERO line, then reaction->delta_g0 is set 
      set to the value on that line, and reaction->deltag0_computed is 
      set to 0. If later compute_reaction_dg0 computes a delta_g0 for the 
      reaction, reaction->deltag0_computed is set to 1.
    */
    reaction->deltag0_computed  = -1;
    forward_rc[0]               = -1.0;
    reverse_rc[0]               = -1.0;
    for (i=0;i<max_regs_per_rxn;i++) {
      reg_constant[reg_base+i] = 1.0;
      reg_exponent[reg_base+i] = 0.0;
      reg_drctn[reg_base+i]    = 1.0;
      reg_species[reg_base+i]  = (int64_t)(-1);
    }
    had_a_title                = 0;
    had_a_dg0                  = 0;
    had_units                  = 0;
    rxn_use                    = 1;
    reg_count                  = 0;
    rxn_ptrs[rxns]              = molecules;
    fgp = fgets(rxn_buffer,rxn_buff_len,rxn_fp);
    line_no = 1;
    state->max_molecule_len = (int64_t)0;
    state->min_molecule_len = rxn_buff_len;
    state->max_compartment_len = (int64_t)0;
    state->min_compartment_len = rxn_buff_len;
    /*
      Build in the empty compartment.
    */
    compartment           = (struct compartment_struct *)&unsorted_cmpts[0];
    compartment->string   = compartment_pos;
    compartment->c_index  = 0;
    compartment->volume   = 0.0;
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
		" %ld characters\n",rxn_buff_len);
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
	  ACTIVITY/ENZYME_LEVEL, PREGULATION and or NREGULATION lines.
	  It also may have a USE_RXN 0 line to turn off a particular
	  reaction for the duration of the run.

      	  Allow pre REACTION lines in reactions file as a header.
      	  So before the first REACTION line the state is -1
      		
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
	  had_a_title = 1;
	  rxn_title_len = (int64_t)(line_len - word1 + 1);
	  /*
	    Copy the reaction title.
	  reaction->title = (char *)&rxn_title_text[rxn_title_pos];
	  */
	  reaction->title = rxn_title_pos;
	  title  = (char *)&rxn_title_text[rxn_title_pos];
	  padding = (align_len - (rxn_title_len & align_mask)) & align_mask;
	  /*
	    Here we should check that rxn_title_pos + 
	    rxn_title_len + padding < max_rxn_title_pos;
	  */
	  rtpd = rxn_title_len + padding;
	  if (rxn_title_pos + rtpd > max_rxn_title_pos) {
	    success = 0;
	    if (lfp) {
	      fprintf(lfp,"parse_reaction_file: Error length of reaction "
		      "titles plus padding > than that measured in "
		      "size_rxns_file.c = %ld\n",
		      max_rxn_title_pos);
	      fflush(lfp);
	    }
	  } else {
	    strncpy(title,(char*)&rxn_buffer[word1],rxn_title_len-1);
	    title[rxn_title_len-1] = '\0';
	  }
	  /*
	    caution bit twiddle follows: 
	    padding = (align_len - (rxn_titl_len % align_len)) % align_len
	    Thou shalt not use %.
	  */
	  rxn_title_pos += rtpd;
	  break;
	case 1:
	  /*
	    Pathway line.
	  */
	  pathway_len = line_len - word1 + 1;
	  padding = (align_len - (pathway_len & align_mask)) & align_mask;
	  /*
	  reaction->pathway = (char *)&pathway_text[pathway_pos];
	  */
	  reaction->pathway = pathway_pos;
	  pathway = (char *)&pathway_text[pathway_pos];
	  ppd = pathway_len + padding;
	  if (pathway_pos + ppd > max_pathway_pos) {
	    success = 0;
	    if (lfp) {
	      fprintf(lfp,"parse_reactions_file: Error pathway text + padding "
		      "length is > than measured by size_rxns_file.c = %ld\n",
		      max_pathway_pos);
	      fflush(lfp);
	    }
	  } else {
	    strncpy(pathway,(char*)&rxn_buffer[word1],word1_len);
	    pathway[word1_len] = '\0';
	  }
	  pathway_pos += (int64_t)ppd;
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
	  reaction->rcompartment = 0;
	  reaction->left_compartment = cmpts;
	  reaction->right_compartment = cmpts;
	  lcompartment = (char *)&compartment_text[compartment_pos];
	  padding = (align_len - ((compartment_len+1) & align_mask)) & align_mask;
	  cpd = compartment_len + 1 + padding;
	  if ((compartment_pos + cpd) > max_compartment_pos) {
	    success = 0;
	    if (lfp) {
	      fprintf(lfp,"parse_reactions_file: Error compartment name "
		      "lengths + padding > that measured in size_rxns_file.c "
		      " = %ld\n",max_compartment_pos);
	      fflush(lfp);
	    }
	  } else {
	    strncpy(lcompartment,(char*)&rxn_buffer[word1],word1_len);
	    lcompartment[word1_len] = '\0';
	    upcase(compartment_len,lcompartment);
	    if ((strcmp(lcompartment,"V") == 0) ||
		(strcmp(lcompartment,"C") == 0)) {
	      if (lfp) {
		fprintf(lfp,"parse_reactions_file: Error line %d: compartments may not have single character names V or C\n",line_no);
		fflush(lfp);
	      }
	      success = 0;
	    }
	  }
	  if (success) {
	    /*
	      unsorted_cmpts->string = lcompartment;
	    */
	    if (cmpts > max_cmpt) {
	      success = 0;
	      if (lfp) {
		fprintf(lfp,"parse_reaction_file: Error parse_reaction parsed"
			" more compartments than size_rxns_file saw, "
			"max_cmpt = %d at line %d\n",max_cmpt,line_no);
		fflush(lfp);
	      }
	    } else {
	      compartment = (struct compartment_struct *)&unsorted_cmpts[cmpts];
	      
	      compartment->string = compartment_pos;
	      compartment->volume = 0.0;
	      compartment->c_index  = cmpts;
	    }
	  }
	  compartment_pos += (int64_t)cpd;
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
	  padding = (align_len - ((compartment_len+1) & align_mask)) & align_mask;
	  cpd = compartment_len + 1 + padding;
	  /*
	  reaction->lcompartment = (char *)&compartment_text[compartment_pos];
	  */
	  reaction->lcompartment = compartment_pos;
	  reaction->left_compartment = cmpts;
	  lcompartment = (char *)&compartment_text[compartment_pos];
	  if ((compartment_pos + cpd) > max_compartment_pos) {
	    success = 0;
	    if (lfp) {
	      fprintf(lfp,"parse_reactions_file: Error compartment name "
		      "lengths + padding > that measured in size_rxns_file.c "
		      " = %ld\n",max_compartment_pos);
	      fflush(lfp);
	    }
	  } else {
	    strncpy(lcompartment,(char*)&rxn_buffer[word1],word1_len);
	    lcompartment[word1_len] = '\0';
	    upcase(compartment_len,lcompartment);
	    if ((strcmp(lcompartment,"V") == 0) ||
		(strcmp(lcompartment,"C") == 0)) {
	      if (lfp) {
		fprintf(lfp,"parse_reactions_file: Error line %d: compartments may not have single character names V or C\n",line_no);
		fflush(lfp);
	      }
	      success = 0;
	    }
	  }
	  /*
	  unsorted_cmpts->string = lcompartment;
	  */
	  if (success) {
	    if (cmpts > max_cmpt) {
	      success = 0;
	      if (lfp) {
		fprintf(lfp,"parse_reaction_file: Error parse_reaction parsed"
			" more compartments than size_rxns_file saw, "
			"max_cmpt = %d at line %d\n",max_cmpt,line_no);
		fflush(lfp);
	      }
	    } else {
	      compartment = (struct compartment_struct *)&unsorted_cmpts[cmpts];
	      compartment->string = compartment_pos;
	      compartment->volume = 0.0;
	      compartment->c_index  = cmpts;
	    }
	  }
	  compartment_pos += (int64_t)cpd;
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
	  padding = (align_len - ((compartment_len+1) & align_mask)) & align_mask;
	  cpd = compartment_len + 1 + padding;
	  /*
	  reaction->rcompartment = (char *)&compartment_text[compartment_pos];
	  */
	  reaction->rcompartment = compartment_pos;
	  reaction->right_compartment = cmpts;
	  rcompartment = (char *)&compartment_text[compartment_pos];
	  if (compartment_pos + cpd > max_compartment_pos) {
	    success = 0;
	    if (lfp) {
	      fprintf(lfp,"parse_reactions_file: Error compartment name "
		      "lengths + padding > that measured in size_rxns_file.c "
		      " = %ld\n",max_compartment_pos);
	      fflush(lfp);
	    }
	  } else {
	    strncpy(rcompartment,(char*)&rxn_buffer[word1],word1_len);
	    rcompartment[word1_len] = '\0';
	    upcase(compartment_len,rcompartment);
	    if ((strcmp(lcompartment,"V") == 0) ||
		(strcmp(lcompartment,"C") == 0)) {
	      if (lfp) {
		fprintf(lfp,"parse_reactions_file: Error line %d: compartments may not have single character names V or C\n",line_no);
		fflush(lfp);
	      }
	      success = 0;
	    }
	  }
	  /*
	  unsorted_cmpts->string = rcompartment;
	  */
	  if (success) {
	    if (cmpts > max_cmpt) {
	      success = 0;
	      if (lfp) {
		fprintf(lfp,"parse_reaction_file: Error parse_reaction parsed"
			" more compartments than size_rxns_file saw, "
			"max_cmpt = %d at line %d\n",max_cmpt,line_no);
		fflush(lfp);
	      }
	    } else {
	      compartment = (struct compartment_struct *)&unsorted_cmpts[cmpts];
	      compartment->string = compartment_pos;
	      compartment->volume = 0.0;
	      compartment->c_index  = cmpts;
	    }
	  }
	  compartment_pos += (int64_t)cpd;
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
	  had_a_dg0 = 1;
	  ns = sscanf ((char*)&rxn_buffer[word1],"%le",
		       &reaction->delta_g0);
	  if (ns < 1) {
	    title  = (char *)&rxn_title_text[reaction->title];
	    if (lfp) {
	      fprintf(lfp,
		      "parse_reactions_file: Error line %d: malformed DGZERO line"
		      " for reaction %s was\n%s\n",
		      line_no,title,rxn_buffer);
	      fflush(lfp);
	    }
	    success = 0;
	    break;
	  }
	  /*
	    Set the reaction deltag0_computed flag to 0 to show that 
	    a dgzero line was supplied for the reaction in the .dat file.
	  */
	  reaction->deltag0_computed = 0;
	  break;
	case 8:
	  /* 
	     A DGZERO-UNITS line
	  */
	  had_units = 1;
	  sl = word1_len;
	  if (sl < 1) {
	    title  = (char *)&rxn_title_text[reaction->title];
	    if (lfp) {
	      fprintf(lfp,
		    "parse_reactions_file: Error line %d: malformed DGZERO-UNITS line,"
		    " for reaction %s, was %s, using KJ/MOL\n",
		      line_no,title,rxn_buffer);
	      fflush(lfp);
	    }
	    reaction->unit_i = 1;
	  } else {
	    upcase(sl,(char*)&rxn_buffer[word1]);
	    if (strncmp((char*)&rxn_buffer[word1],"KJ/MOL",6) == 0) {
	      reaction->unit_i = 1;
	    } else {
	      reaction->unit_i = 0;
	    }
	  }
	  break;
        case 9:
	case 12:
	  /*
	    An ACTIVITY/ENZYME_LEVEL line.
	  */
	  ns = sscanf ((char*)&rxn_buffer[word1],"%le",
		       &reaction->enzyme_level);
	  if (ns < 1) {
	    if (lfp) {
	      fprintf(lfp,
		      "parse_reactions_file: Error line %d: malformed "
		      "ACTIVITY/ENZYME_LEVEL line was\n%s\n",
		      line_no,rxn_buffer);
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
	      fprintf(lfp,"parse_reactions_file: Error line %d: limit of %d "
		      "regulations per reaction was exceeded. Increase the "
		      "MAX_REGS_PER_RXN parameter and rerun.\n",
		      line_no,max_regs_per_rxn);
	      fflush(lfp);
	    }
	    break;
	  }
	  metabolite = (char *)&regulation_text[regulation_pos];
	  padding = (align_len - ((word1_len+1) & align_mask)) & align_mask;

	  rgpd = word1_len + 1 + padding;
	  if ((regulation_pos + rgpd) > max_regulation_pos) {
	    success = 0;
	    if (lfp) {
	      fprintf(lfp,"parse_reactions_file: Error regulation text "
		      "lengths + padding > that measured in size_rxns_file.c "
		      " = %ld\n",max_regulation_pos);
	      fflush(lfp);
	    }
	  } else {
	    strncpy(metabolite,(char*)&rxn_buffer[word1],word1_len);
	    metabolite[word1_len] = '\0';
	    upcase(word1_len,metabolite);
	    /*
	      reg_pos = reg_base + reg_count;
	    */
	    reg_species[reg_pos] = regulation_pos;
	    reg_drctn[reg_pos]   = 1.0;
	    nr = sscanf((char*)&rxn_buffer[word1+word1_len],"%le %le",&reg_constant[reg_pos],&reg_exponent[reg_pos]);
	    if (nr < 2) {
	      success = 0;
	      if (lfp) {
		fprintf(lfp,"parse_reactions_file: Error line %d: PREGULATION LINE, %s\n"
		      "did not have valid constant and exponent.\n",
			line_no,rxn_buffer);
		fflush(lfp);
	      }
	    }
	  }
	  regulation_pos += rgpd;
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
	      fprintf(lfp,"parse_reactions_file: Error line %d: limit of %d "
		      "regulations per reaction was exceeded. Increase the "
		      "MAX_REGS_PER_RXN parameter and rerun.\n",
		      line_no,max_regs_per_rxn);
	      fflush(lfp);
	    }
	    break;
	  }
	  metabolite = (char *)&regulation_text[regulation_pos];
	  padding = (align_len - ((word1_len+1) & align_mask)) & align_mask;

	  rgpd = word1_len + 1 + padding;
	  if ((regulation_pos + rgpd) > max_regulation_pos) {
	    success = 0;
	    if (lfp) {
	      fprintf(lfp,"parse_reactions_file: Error regulation text "
		      "lengths + padding > that measured in size_rxns_file.c "
		      " = %ld\n",max_regulation_pos);
	      fflush(lfp);
	    }
	  } else {
	    strncpy(metabolite,(char*)&rxn_buffer[word1],word1_len);
	    metabolite[word1_len] = '\0';
	    upcase(word1_len,metabolite);
	    /*
	      reg_pos = reg_base + reg_count;
	    */
	    reg_species[reg_pos] = regulation_pos;
	    reg_drctn[reg_pos]   = 0.0;
	    nr = sscanf((char*)&rxn_buffer[word1+word1_len],"%le %le",&reg_constant[reg_pos],&reg_exponent[reg_pos]);
	    if (nr < 2) {
	      success = 0;
	      if (lfp) {
		fprintf(lfp,"parse_reactions_file: Error line %d: NREGULATION "
			"LINE, %s\n did not have valid constant and exponent.\n",
			line_no,rxn_buffer);
		fflush(lfp);
	      }
	    }
	  }
	  regulation_pos += rgpd;
	  reg_count += 1;
	  reg_pos += 1;		 
	  break;
	case 13:
	  /*
	    // reaction terminator line
	  */
	  /*
	    First we should check to see that a left
	    and a right side and a title line were seen.
	  */
	  num_reactants = reaction->num_reactants;
	  num_products  = reaction->num_products;
	  if ((num_reactants == 0) || (num_products == 0) ||
	      (had_a_title == 0)) {
	    success = 0;
	    /*
	      Incomplete reaction specification ignore it.
	    */
	    if (lfp) {
	      fprintf(lfp,"parse_reactions_file: Error: Incomplete reaction on"
		      " reaction ending on line %d.\n"
		      "Each reaction must have REACTION, LEFT, RIGHT, DGZER0,"
		      " and DGZERO-UNITS lines\n",line_no);
	      fflush(lfp);
	    }
	    break;
	  } else {
	    reaction->self_id = rxns;
	    activities[rxns]  = reaction->enzyme_level;
	    enzyme_level[rxns] = reaction->enzyme_level;
	    forward_rc[rxns]   = reaction->forward_rc;
	    reverse_rc[rxns]   = reaction->reverse_rc;
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
	    sum_coeffs = 0;
	    for (j = mol_pos;j<mol_pos_lim;j++) {
	      sum_coeffs += coefficients[j];
	      /*
		Only look to set the compartment index if the molecule did
		not have a local compartment (:compartment).
	      */
	      if (rxn_molecules->c_index == -1) {
		if (coefficients[j] < 0) {
		  rxn_molecules->c_index = reaction->left_compartment;
		} else {
		  rxn_molecules->c_index = reaction->right_compartment;
		}
	      }
	      rxn_molecules += 1; /* Caution address arithmetic */
	    }
	    reaction->coefficient_sum = sum_coeffs;
	    coeff_sum[rxns] = sum_coeffs;
	    use_rxn[rxns]   = rxn_use;
	    rxns += 1;
	    /*
	      Caution address arithmetic follows
	      reaction = (struct reaction_struct*)&reactions[rxns];
	    */
	    reaction += 1;
	    reg_base  += max_regs_per_rxn;
	  }
	  reg_pos   =  reg_base;
	  reg_count = 0;
	  had_a_title = 0;
	  had_a_dg0   = 0;
	  had_units   = 0;
	  if (rxns < (int)state->number_reactions) {
	    reaction->lcompartment      = 0;
	    reaction->rcompartment      = 0;
	    reaction->pathway           = -1;
	    reaction->left_compartment  = 0;
	    reaction->right_compartment = 0;
	    reaction->num_reactants     = 0;
	    reaction->num_products      = 0;
	    reaction->activity          = 1.0;
	    reaction->enzyme_level      = 1.0;
	    reaction->forward_rc        = -1.0;
	    reaction->reverse_rc        = -1.0;
	    reaction->unit_i            = 1;
	    reaction->delta_g0          = 0.0;
	    forward_rc[rxns]            = -1.0;
	    reverse_rc[rxns]            = -1.0;
	    /*
	      The following three lines added by DGT on 4/18/2013
	    */
	    reaction->temp_kelvin       = state->temp_kelvin;
	    reaction->ph                = state->ph;
	    reaction->ionic_strength    = state->ionic_strength;

	    reaction->deltag0_computed  = -1;
	    
	    rxn_ptrs[rxns]              = molecules;
	    rxn_use                     = 1;
	    for (i=0;i<max_regs_per_rxn;i++) {
	      reg_constant[reg_base+i] = 1.0;
	      reg_exponent[reg_base+i] = 0.0;
	      reg_drctn[reg_base+i]    = 1.0;
	      reg_species[reg_base+i]  = (int64_t)(-1);
	    }
	  }
	  break;
	case 14: /* K_FORWARD */
	  /*
	    A FOWARD rate line.
	  */
	  ns = sscanf ((char*)&rxn_buffer[word1],"%le",
		       &reaction->forward_rc);
	  if (ns < 1) {
	    if (lfp) {
	      fprintf(lfp,
		     "parse_reactions_file: Error line %d: malformed "
		      "FORWARD_RATE line was\n%s\n",
		      line_no,rxn_buffer);
	      fflush(lfp);
	    }
	    success = 0;
	    break;
	  }
	  break;
	case 15: /* K_REVERSE */
	  /*
	    A REVERSE rate line.
	  */
	  ns = sscanf ((char*)&rxn_buffer[word1],"%le",
		       &reaction->reverse_rc);
	  if (ns < 1) {
	    if (lfp) {
	      fprintf(lfp,
		      "parse_reactions_file: Error line %d: malformed "
		      "REVERSE_RATE line was\n%s\n",
		      line_no,rxn_buffer);
	      fflush(lfp);
	    }
	    success = 0;
	    break;
	  }
	  break;
	case 16:
	  /*
	    A USE_RXN line.
	  */
	  ns = sscanf((char*)&rxn_buffer[word1],"%d",&rxn_use);
	  if (rxn_use == 0) {
	    rxn_not_used = 1;
	  }
	  break;
	case 17:
	  /*
	    An IGNORE line.
	  */
	  rxn_use = 0;
	  rxn_not_used = 1;
	  break;
        default:
	break;
	} /* end switch(line_type) */
      }/* end if (success) */
      fgp = fgets(rxn_buffer,rxn_buff_len,rxn_fp);
      line_no += 1;
    } /* end while(fgp...) */
    rxn_ptrs[rxns] = molecules;
    /*
      Check that last line was a //.
    */
    /* 
       this comes from the postion of //  in the keywords list in 
       init_rxn_file_keywords.c 
    */
    last_line_type = 13; 
    if (success) {
      if (line_type != last_line_type) {
	if (lfp) {
	  fprintf(lfp,
		  "parse_reactions_file: Error line %d: reactions file did "
		  "not end in //\n",line_no);
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
  /*
    If one or more reactions was turned off with a USE_RXN 0 line, or an 
    ignore line, we need to remove that reaction from the reactions struct,
    and from the reactions matrix and adjust the num_reactions
    field.
  */
  if (success) {
    if (rxn_not_used) {
      success = boltzmann_compress_reactions(state);
    }
  }
  if (success) {
    state->number_compartments = cmpts;
    state->number_molecules    = molecules;
  }
  return(success);
}
