/* rxn_map_run.c
*******************************************************************************
rxn_mapp

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


#ifdef TIMING_ON
struct timing_struct timing_data;
#endif

#ifdef LIBUNWIND
#include "luwtb.h"
#endif
/*
#define RXN_MAP_DBG 1
*/
#include "rxn_map_run.h"
int rxn_map_run(struct state_struct *state,
		int mi, int mj) {
  /*
    Called by: rxn_map
  */
  struct molecules_matrix_struct *molecules_matrix;
  struct reactions_matrix_struct *rxn_matrix;
  struct molecule_struct *sorted_molecules;
  struct compartment_struct *sorted_compartments;
  struct molecule_struct *molecule_i;
  struct molecule_struct *molecule_j;
  struct compartment_struct *compartment_i;
  struct compartment_struct *compartment_j;
  struct stack_level_elem_struct sles;
  struct stack_level_elem_struct *levels;
  struct stack_level_elem_struct *cur_level;
  struct stack_level_elem_struct *next_level;
  struct stack_level_elem_struct *level_i;
  char   *molecules_text;
  char   *compartment_text;

  int64_t *rxn_ptrs;
  int64_t *m_indices;
  int64_t *r_coefs;

  int64_t *mol_ptrs;
  int64_t *r_indices;
  int64_t *m_coefs;


  int64_t *r_unmarked;
  int64_t *stack;
  int64_t *terminal_reaction;
  int64_t *terminal_direction;
  int64_t stack_pos_lim;
  int64_t raw_path_count;
  int64_t stack_pos;
  int64_t stack_entry_size;
  int64_t rxn_l;
  int64_t rxn_i;
  int64_t mol_j;
  int64_t rxn_k;
  int64_t rxn_j;
  int64_t dxn_i;
  int64_t rxn_cand;
  
  int64_t num_reactions;
  int64_t nunique_molecules;
  int64_t k;
  int64_t i;
  int64_t j;
  int64_t ci;
  int64_t cj;
  int64_t ask_for;
  int64_t one_l;
  int64_t zero_l;
  int64_t nzr;
  int64_t direction;
  int success;
  int level;
  /*
    Calls:
  */
  /*
    so we want to generate reaction chains that start with one of the
    fixed concentration products, and end with another one of the fixed
    concentration products, with no cycles and each rxn used in only
    one direction in the chain. Then we will want to form the summary 
    reaction for each chain, printing that out. In theory then
    any molecule appearing in the summary reaction for some chain should
    also be a fixed concentration molecule. So ultimately we want to
    produce a list of the fixed concentration molecules.

    Getting started, 


  */
  molecules_matrix    = state->molecules_matrix;
  rxn_matrix          = state->reactions_matrix;
  sorted_molecules    = state->sorted_molecules;
  sorted_compartments = state->sorted_compartments;
  nunique_molecules   = state->nunique_molecules;
  num_reactions       = state->number_reactions;
  nzr                 = state->number_molecules;
  molecules_text      = state->molecules_text;
  compartment_text    = state->compartment_text;

  rxn_ptrs            = rxn_matrix->rxn_ptrs;
  m_indices           = rxn_matrix->molecules_indices;
  r_coefs             = rxn_matrix->coefficients;

  mol_ptrs            = molecules_matrix->molecules_ptrs;
  r_indices           = molecules_matrix->reaction_indices;
  m_coefs             = molecules_matrix->coefficients;

  zero_l = (int64_t)0;
  one_l  = (int64_t)1;
  stack_pos_lim       = (num_reactions * num_reactions) + num_reactions;
  ask_for = stack_pos_lim * ((int64_t)sizeof(int64_t));
  success = 1;
  raw_path_count = 0;
  stack = (int64_t*)calloc(one_l,ask_for);
  if (stack == NULL) {
    fprintf(stderr,"rxn_map_run: Could not allocate %ld bytes for stack\n",
	    ask_for);
    fflush(stderr);
    success = 0;
  }
  if (success) {
    ask_for = num_reactions * ((int64_t)sizeof(sles));
    levels = (struct stack_level_elem_struct *)calloc(one_l,ask_for);
    if (levels == NULL) {
      fprintf(stderr,"rxn_map_run: Could not allocate %ld bytes for levels\n",
	      ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  if (success) {
    ask_for = num_reactions * ((int64_t)sizeof(int64_t));
    r_unmarked = (int64_t*)calloc(one_l,ask_for);
    if (r_unmarked == NULL) {
      fprintf(stderr,
	      "rxn_map_run: Could not allocate %ld bytes for r_unmarked\n",
	      ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  if (success) {
    ask_for = num_reactions * ((int64_t)sizeof(int64_t));
    terminal_reaction = (int64_t*)calloc(one_l,ask_for);
    if (terminal_reaction == NULL) {
      fprintf(stderr,
      "rxn_map_run: Could not allocate %ld bytes for terminal_reaction\n",
	      ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  if (success) {
    ask_for = num_reactions * ((int64_t)sizeof(int64_t));
    terminal_direction = (int64_t*)calloc(one_l,ask_for);
    if (terminal_direction == NULL) {
      fprintf(stderr,
      "rxn_map_run: Could not allocate %ld bytes for terminal_direction\n",
	      ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  /*
  direction_list      = &path_list[]
  */
  
  if (success) {
    molecule_i = (struct molecule_struct *)&sorted_molecules[mi];
    ci = molecule_i->c_index;
    compartment_i = (struct compartment_struct *)&sorted_compartments[ci];
    molecule_j = (struct molecule_struct *)&sorted_molecules[mj];
    cj = molecule_j->c_index;
    compartment_j = (struct compartment_struct *)&sorted_compartments[cj];
    /* 
      Now we want to build chains of reactions, that
      start with mi and end with mj.
      Total length of the chains must be less than
      num_reactions.
      Notion of a current level, unprocessed reactions at this level,
      Each level has a current reaction index, a last reaction index
      into the stack (int array).
      current reaction being processed at this level,
      look at his "product molecules" walk through their reaction lists,
      adding reactions not yet in this chain and not yet marked for this level.

      Random thoughts, one might look at reaction pairs,
      reactions that have a common molecule on opposite sides of their
      equation, then use a dynamic programming approach to stitch
      these together matching ends that are the same, looking at all
      2 length reactions, then all 3 length reactions, (merging two length
      reactions with a matching right/left end), then 4 length reactions,
      merging a 3 with a 2 or a 2 with a 3, building successively longer
      chains
      
    */
    fprintf (stdout,"mi = %d, mj = %d\n",mi,mj);
    if (ci > 0) {
      fprintf(stdout,"%s:%s to ",
	      (char*)&molecules_text[molecule_i->string],
	      (char*)&compartment_text[compartment_i->string]);
    } else {
      fprintf(stdout,"%s to ",
	      (char*)&molecules_text[molecule_i->string]);
    }
    if (cj > 0) {
      fprintf(stdout,"%s:%s\n",
	      (char*)&molecules_text[molecule_j->string],
	      (char*)&compartment_text[compartment_j->string]);
    } else {
      fprintf(stdout,"%s\n",
	      (char*)&molecules_text[molecule_j->string]);
    }
    fflush(stdout);
    /*
      We will stack reactions to be examined on layers, (distance
      from initial reaction. Each layer (max num layers = num_reactions)
      will have a current index ptr, and a last_ptr, when we have exhausted
      the last reaction on a level that level is "popped" and we advance
      to the next reaction in the layer below. It will be useful to
      have a indicator array for possible terminal reactions so we 
      no when to record a complete path. Also each reaction needs to
      have a direction associated with it, because we start at rxn 0
      we need to keep a separate direction indicator (we could start
      at one though but that has complications, we'll think about that).
      So the stack size will need to be num_reactions*num_reactions + 
      num_reactions long as 
      each level can have at most num_reactions-level additional paths
      to investigate, hence num_reactions * (num_reactions+1)/2 entries
      with each entry having a reaction number and a direction to store.
    */
    /*
      Mark the terminal reactions.
    */
    for (i=0;i<num_reactions;i++) {
      terminal_reaction[i] = zero_l;
    }
    for (j=mol_ptrs[mj];j<mol_ptrs[mj+1];j++) {
      rxn_j = r_indices[j];
      terminal_reaction[rxn_j] = one_l;
      terminal_direction[rxn_j] = one_l;
      if (m_coefs[j] < 0) {
	terminal_direction[rxn_j] = -one_l;
      }
    }
    stack_pos  = zero_l;
    /*
      Unmark all the reactions.
    */
    for (k=0;k<num_reactions;k++) {
      r_unmarked[k] = one_l;
    }
    /*
      Add the initial reactions to the stack.
    */
    stack_entry_size = one_l + one_l;
    for (i=mol_ptrs[mi];i<mol_ptrs[mi+1];i++) {
      rxn_cand = r_indices[i];
      if (r_unmarked[rxn_cand]) {
	r_unmarked[rxn_cand] = 0;
	stack[stack_pos] = (int)rxn_cand;
	stack[stack_pos+1] = one_l;
	if (m_coefs[i] < (int64_t)0) {
	  stack[stack_pos+1] = - one_l;
	}
	stack_pos += stack_entry_size;
      }
    }
    /*
      Unmark all the reactions.
    */
    for (k=0;k<num_reactions;k++) {
      r_unmarked[k] = 1;
    }
    level = 0;
    cur_level = (struct stack_level_elem_struct*)&levels[level];
    cur_level->first_stack_pos  = 0;
    cur_level->last_stack_pos   = (mol_ptrs[mi+1] - mol_ptrs[mi] - 1) << 1;
    cur_level->cur_stack_pos    = cur_level->first_stack_pos;    
    next_level = (struct stack_level_elem_struct*)&levels[level+1];
    next_level->first_stack_pos = cur_level->last_stack_pos + stack_entry_size;
    stack_pos = cur_level->first_stack_pos;
    raw_path_count      = 0;
    while (level >= 0) {
      cur_level = (struct stack_level_elem_struct*)&levels[level];
      /*
	Check to see if we are finished with this
	level.
      */
      if ((cur_level->cur_stack_pos) > (cur_level->last_stack_pos)) {
	/*
	  decrement level 
	  advance level_cur for previous level.
	*/
	level -= 1;
	if (level >= 0) {
	  cur_level = (struct stack_level_elem_struct*)&levels[level];
	  stack_pos = cur_level->cur_stack_pos;
	  rxn_l = stack[stack_pos];
	  r_unmarked[rxn_l] = 1;
	  cur_level->cur_stack_pos += stack_entry_size;
	}
      } else {
	stack_pos = cur_level->cur_stack_pos;
	rxn_l     = stack[stack_pos];
	direction = stack[stack_pos+1];
	next_level = (struct stack_level_elem_struct*)&levels[level+1];
	next_level->first_stack_pos = cur_level->last_stack_pos + 
	                              stack_entry_size;
	r_unmarked[rxn_l] = 0;
	/*
	  Check to see if this reaction is a terminal reaction,
	  If so print out the reaction chain.
	*/
	if ((terminal_reaction[rxn_l]) && 
	    (terminal_direction[rxn_l] == direction)) {
	  raw_path_count += 1;
	  for (i=0;i<level;i++) {
	    level_i = (struct stack_level_elem_struct*)&levels[i]; 
	    stack_pos = level_i->cur_stack_pos;
	    rxn_i = stack[stack_pos];
	    dxn_i = stack[stack_pos+1];
	    if (dxn_i < 0) {
	      fprintf(stdout,"-%ld\t",rxn_i);
	    } else {
	      fprintf(stdout,"%ld\t",rxn_i);
	    }
	  }
	  if (direction < 0) {
	    fprintf(stdout,"-%ld\n",rxn_l);
	  } else {
	    fprintf(stdout,"%ld\n",rxn_l);
	  }
	}
	/*
	  Mark the reactions from levels in the chain.
	  This might not be necesarry.
	for (i=0;i<level;i++) {
	  stack_pos = level_cur[level];
	  rxn_i     = stack[stack_pos];
	  r_unmarked[rxn_i] = 0;
	}
	*/
	/*
	  Walk through molecules on the production side of rxn_l.
	*/
	if (level < num_reactions - 1) {
	  for (j=rxn_ptrs[rxn_l];j<rxn_ptrs[rxn_l+1];j++) {
	    mol_j = m_indices[j];
	    if ((r_coefs[j] * direction) > 0) {
	      /*
	        Molecule is produced by this reaction in that direction.
	        so walk through its list of reactions adding those that 
	        are unmarked to the stack.
	      */
	      level +=1;
	      stack_pos = next_level->first_stack_pos;
	      for (k=mol_ptrs[mol_j];k<mol_ptrs[mol_j+1];k++) {
	        rxn_k = r_indices[k];
	        if (r_unmarked[rxn_k]) {
		  r_unmarked[rxn_k] = 0;
		  stack[stack_pos] = rxn_k;
		  stack[stack_pos+1] = 1;
		  if (m_coefs[k] > 0) {
		    stack[stack_pos+1] = -1;
		  }
		  stack_pos += stack_entry_size;
	        }
	      }
	      next_level->last_stack_pos = stack_pos - stack_entry_size;
	      next_level->cur_stack_pos  = next_level->first_stack_pos;
	      /*
	        Now want to unmark the reactions just added to the stack.
	      */
	      for (k=next_level->first_stack_pos;
		   k<=next_level->last_stack_pos;
		   k+=stack_entry_size) {
	        r_unmarked[stack[k]] = 1;
	      }
	    } /* end if (molecule was a product) */
	  } /* end for (j...) loop over molecules */
	} else {
	  /* last level */
	  cur_level->cur_stack_pos += stack_entry_size;
	}
      } /* end else not the end of a level */
    } /* end while level */
  } /* end if allocate ws succeeded */
  fprintf(stdout,"raw_path_count  = %ld\n",raw_path_count);
  fflush(stdout);
  return (success);
}
	 
