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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <float.h>
#include <signal.h>
#include <unistd.h>

#include "djb_timing_b.h"
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
		struct molecules_matrix_struct *molecules_matrix,
		int i, int j) {
  struct rxn_matrix_struct *rxn_matrix;
  struct istring_elem_struct *sorted_molecules;
  struct istring_elem_struct *sorted_compartments;
  struct istring_elem_struct *molecule_i;
  struct istring_elem_struct *molecule_j;
  struct istring_elem_struct *molecule_t;
  struct istring_elem_struct *compartment_i;
  struct istring_elem_struct *compartment_j;
  char   *molecules_text;
  char   *compartment_text;

  int64_t *ws;
  int64_t *rxn_ptrs;
  int64_t *m_indices;
  int64_t *r_coefs;

  int64_t *mol_ptrs;
  int64_t *r_indices;
  int64_t *m_coefs;


  int64_t *m_unmarked;
  int64_t *r_unmarked;
  int64_t *stack;
  int64_t *path_ptrs;
  int64_t *path_list;
  int64_t stack_pos_lim;
  int64_t path_count_lim;
  int64_t raw_path_count;
  int64_t path_pos;
  int64_t stack_pos;
  
  int64_t num_reactions;
  int64_t nunique_molecules;
  int64_t k;
  int64_t mi;
  int64_t mj;
  int64_t mt;
  int64_t ci;
  int64_t cj;
  int64_t ask_for;
  int64_t one_l;
  int64_t nzr;
  int64_t path_pos_lim;
  int64_t c_molecule;
  int64_t next_reaction_pos;
  int64_t potential_reaction;
  int64_t direction;
  int64_t potential_molecule;
  int64_t next_molecule_pos;
  int success;
  int padi;
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
  rxn_matrix          = state->reactions_matrix;
  sorted_molecules    = state->sorted_molecules;
  sorted_compartments = state->sorted_cmpts;
  nunique_molecules   = state->nunique_molecules;
  num_reactions       = state->number_reactions;
  nzr                 = state->number_molecules;
  molecules_text      = state->molecules_text;
  compartment_text    = state->compartment_text;

  rxn_ptrs            = rxn_matrix->rxn_ptrs;
  m_indices           = rxn_matrix->molecules_indices;
  r_coefs             = rxn_matrix->coefficients;

  mol_ptrs            = molecules_matrix->molecules_ptrs;
  r_indices           = molecules_matrix->rxn_indices;
  m_coefs             = molecules_matrix->coefficients;

  /*
    Comment here then ws needs to be of size 
    (nunique_molecules + num_reactions + stack_pos_lim + path_count_lim
    + path_pos_lim.
  */
  one_l = (int64_t)1;
  stack_pos_lim       = (num_reactions << 3) + 4;
  path_count_lim      = nzr*nzr*nzr;
  path_pos_lim        = path_count_lim * num_reactions;
  ask_for = (nunique_molecules + num_reactions + stack_pos_lim + path_count_lim
	     + path_count_lim + path_pos_lim);
  ask_for = ask_for*sizeof(int64_t);
  ws = (int64_t*)calloc(one_l,ask_for);
  if (ws == NULL) {
    fprintf(stderr,"rxn_map_run: Could not allocate %ld bytes for ws\n",
	    ask_for);
    fflush(stderr);
    success = 0;
  }
  m_unmarked          = &ws[0];
  r_unmarked          = &ws[nunique_molecules];
  stack               = &r_unmarked[num_reactions];
  raw_path_count      = 0;
  path_pos            = 0;
  path_ptrs           = &stack[stack_pos_lim];
  path_list           = &path_ptrs[path_count_lim];
  direction_list      = &path_list[
  
  if (success) {
    /*
      So we need a loop over the sorted molecules.
    */
    path_ptrs[0] = 0;
    molecule_i = (struct istring_elem_struct *)&sorted_molecules[i];
    mi = molecule_i->m_index;
    ci = molecule_i->c_index;
    compartment_i = (struct istring_elem_struct *)&sorted_compartments[ci];
    molecule_j = (struct istring_elem_struct *)&sorted_molecules[j];
    mj = molecule_j->m_index;
    cj = molecule_j->c_index;
    compartment_j = (struct istring_elem_struct *)&sorted_compartments[cj];
    /* 
      Now we want to build chains of reactions, that
      start with mi and end with mj where mj > mi and mj is fixed.
      Total length of the chains must be less than
      nunique_molecules^2.
    */
    fprintf (stdout,"mi = %ld, mj = %ld\n",mi);
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
    stack_pos  = 0;
    stack[0] = mi;
    stack[1] = 0;
    stack[2] = mol_ptrs[mi];
    stack[3] = mol_ptrs[mi+1];
    /* Top level start molecule. Unmark all other molecules 
       and all reactions.
    */
    for (k=0;k<num_reactions;k++) {
      r_unmarked[k] = 1;
    }
    for (k=0;k<nunique_molecules;k++) {
      m_unmarked[k] = 1;
    }
    m_unmarked[mi] = 0;
    while (stack_pos >= 0) { 
      if ((stack_pos & 7) == 0) {
	/*
	  A molecule is on top of the stack.
	  First check to see if it is a desired end molecule mj.
	*/
	mt = stack[stack_pos];
	molecule_t = (struct istring_elem_struct *)&sorted_molecules[mt];
	
	if (mt == mj) {
	  for (k=4;k<stack_pos;k+=8) {
	    if (path_pos >= path_pos_lim) {
	      fprintf(stderr,
		      "rxn_map_run: path_pos_lim of %ld reached\n",
		      path_pos_lim);
	      success = 0;
	      stack_pos = -4;
	      break;
	    } else {
	      path_list[path_pos] = stack[k];
	      direction_list[path_pos] = stack[k+1];
	      if (stack[k+1] < 0) {
		fprintf(stdout,"-%ld\n",path_list[path_pos]);
	      } else {
		fprintf(stdout,"%ld\n",path_list[path_pos]);
	      }
	      path_pos += 1;
	    }
	  } /* end for(k...) */
	  fprintf(stdout,"-----------------\n");
	  fflush(stdout);
	  raw_path_count += 1;
	  if (raw_path_count >= path_count_lim) {
	    fprintf(stderr,"rxn_map_run: raw_path_count exceeded %ld\n",
		    path_count_lim);
	    fflush(stderr);
	    success = 0;
	  }
	  path_ptrs[raw_path_count] = path_pos;
	}
	if (success == 0) {
	  break;
	}
	next_reaction_pos = stack[stack_pos+2];
	if (next_reaction_pos >= stack[stack_pos+3]) {
	  /* 
	    No more reactions to process for this molecule,
	    pop the molecule, and unmark it.
	  */
	  if (stack_pos == 0) {
	    break;
	  }
	  m_unmarked[stack[stack_pos]] = 1;
	  stack_pos -= 4;
	} else {
	  stack[stack_pos+2] += 1;
	  potential_reaction = r_indices[next_reaction_pos];
	  if (r_unmarked[potential_reaction]) {
	    r_unmarked[potential_reaction] = 0;
	    if (m_coefs[next_reaction_pos] < 0) {
	      direction = (int64_t)1;
	    } else {
	      direction = (int64_t)(-1);
	    }
	    stack_pos += 4;
	    if (stack_pos > stack_pos_lim) {
	      fprintf(stderr,"stack_pos exceeded %ld\n",
		      stack_pos_lim);
	      fflush(stderr);
	      stack_pos = -4;
	      success = 0;
	      break;
	    } else {
	      stack[stack_pos]   = potential_reaction;
	      stack[stack_pos+1] = direction;
	      stack[stack_pos+2] = rxn_ptrs[potential_reaction];
	      stack[stack_pos+3] = rxn_ptrs[potential_reaction+1];
	    }
	  } /* end if unmarked reaction */
	} /* end else pushed a reaction. */
      } else {
	/*
	  stack_pos is not a multiple of eight, therefore it
	  points to a reaction node.
	  Look for an unmarked molecule on the 
	  production side of the reaction to add to the chain.
	*/
	next_molecule_pos = stack[stack_pos+2];
	if (next_molecule_pos >= stack[stack_pos+3]) {
	  r_unmarked[stack[stack_pos]] = 1;
	  stack_pos -= 4;
	} else {
	  stack[stack_pos+2] += 1;
	  direction = stack[stack_pos+1];
	  potential_molecule = m_indices[next_molecule_pos];
	  /* 
	    Check for an unmarked molecule.
	  */
	  if (m_unmarked[potential_molecule]) {
	    if (r_coefs[next_molecule_pos] * direction > 0) {
	      m_unmarked[potential_molecule] = 0;
	      stack_pos += 4;
	      if (stack_pos > stack_pos_lim) {
		fprintf(stderr,"stack_pos exceeded %ld\n",
			stack_pos_lim);
		fflush(stderr);
		stack_pos = -4;
		success = 0;
		break;
	      } else {
		stack[stack_pos]   = potential_molecule;
		stack[stack_pos+1] = 0;
		stack[stack_pos+2] = mol_ptrs[potential_molecule];
		stack[stack_pos+3] = mol_ptrs[potential_molecule+1];
	      }
	    }/* end unmarked molecule on product side */
	  } /* end unmarked_molecule */
	} /* end not the end of the molecule list. */
      } /* end else a reaction node */
    } /* end while stack_pos >= 0 */
  } /* end if allocate ws succeeded */
  fprintf(stdout,"raw_path_count  = %ld\n",raw_path_count);
  fflush(stdout);
  return (success);
}
	 
