/* parse_side_line.c
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

#include "count_ws.h"
#include "count_nws.h"
#include "upcase.h"
#include "is_a_coef.h"
#include "find_colon.h"

#include "parse_side_line.h"
int parse_side_line(char *rctnts_p, 
		    int64_t *molecules_pos_p,
		    int64_t *compartment_pos_p,
		    int *molecules_p,
		    int *cmpts_p,
		    struct rxn_struct *reaction,
		    struct state_struct *state,
		    int side) {
  /*
    Parse a left/right line, count and record reactant molecules and coefficients.
    Variable           TMF    Description

    rctnts_p           C*I   Pointer to the string of reactants/products.
    molecules_pos_p    L*B   Index into the raw_molecules_text and 
                             molecules_text buffers in the state struct.
			     The contents of this pointer is advanced 
			     in this routine as molecules are encountered.
		       	   
    compartment_pos_p  L*B   Index into the compartment_text buffer
                             in the state struct. The contents of
			     this is pointer is advanced in this routine as 
			     compartments are encountered.

    molecules_p        I*B   Address of the count of molecules encountered
                             so far. Contents are updated by this routine.

    cmpts_p            I*B   Address of the count of compartments encountered
                             so far. Contents are updated by this routine.
 
    reaction           G*I   Pointer to the current reaction data struct.
                             The num_reactants field or num_products field
			     of this struct is incremented depending on whether
			     side is -1 or 1.

    state              G*I   Pointer to the global state struct.
                             The unsorted_molecules and 
			     unsorted_compartments arrays are modified.
			     The compartment_text, raw_molecules_text,
			     and molecules_text fields are also updated.
			     align_mask and align_len fields are inputs
			     but not modified.
			     The reactions_matrix struct is updated,
			     specifically its molecules_indices, coefficients,
			     and matrix_text fields.
			     
    side               ISI   -1 for a left line, 1 for a right line.

    
  */
  struct istring_elem_struct *unsorted_molecules;
  struct istring_elem_struct *unsorted_cmpts;
  struct rxn_matrix_struct *rxns_matrix;
  int64_t *molecules_indices;
  int64_t *coefficients;
  int64_t *matrix_text;
  int64_t molecules_pos;
  int64_t compartment_pos;
  int64_t align_len;
  int64_t align_mask;
  int64_t len;
  int64_t pos;
  int64_t sll;

  char    *rctnts;
  char    *compartment_text;
  char    *molecules_text;
  char    *raw_molecules_text;
  char    *compartment;

  int     sl;
  int     colon_loc;

  int     skip;
  int     ci;

  int     padding;
  int     compartment_len;

  int     cmpts;
  int     molecules;

  int     success;
  int     padi;
  success            = 1;
  molecules_pos      = *molecules_pos_p;
  compartment_pos    = *compartment_pos_p;
  molecules          = *molecules_p;
  cmpts              = *cmpts_p;
  unsorted_molecules = (struct istring_elem_struct *)&state->unsorted_molecules[molecules];
  unsorted_cmpts     = (struct istring_elem_struct *)&state->unsorted_cmpts[cmpts];
  rctnts             = rctnts_p;

  rxns_matrix        = state->reactions_matrix;
  molecules_indices  = rxns_matrix->molecules_indices;
  coefficients       = rxns_matrix->coefficients;
  matrix_text        = rxns_matrix->text;

  align_len    	     = state->align_len;
  align_mask   	     = state->align_mask;
  compartment_text   = state->compartment_text;
  molecules_text     = state->molecules_text;
  raw_molecules_text = state->raw_molecules_text;

  pos = 0;
  len = strlen(rctnts);
  while (pos < len) {
    sl = count_nws(rctnts);
    if (sl > 0) {
      molecules_indices[molecules] = molecules;
      if (is_a_coef(sl,rctnts)) {
	rctnts[sl] = '\0';
	if (side < 0) {
	  coefficients[molecules] = - atoi(rctnts);
	} else {
	  coefficients[molecules] = atoi(rctnts);
	}
	rctnts[sl] = ' ';

	pos += (int64_t)sl;
	rctnts += sl; /* Caution address arithmetic. */

	skip   =  count_ws(rctnts);

	pos += skip;
	rctnts += skip; /* Caution address arithmetic. */

	sl = count_nws(rctnts);
      } else {
	coefficients[molecules] = (double)side;
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
	/*
	unsorted_cmpts->string = (char *)&compartment_text[compartment_pos];
	*/
	unsorted_cmpts->string = compartment_pos;
	compartment = (char *)&compartment_text[compartment_pos];
	unsorted_cmpts->c_index  = cmpts;
	unsorted_cmpts += 1; /* Caution address arithmetic */
	ci = cmpts;
	cmpts += 1;
	strcpy(compartment,(char*)&rctnts[colon_loc+1]);
	upcase(compartment_len,compartment,
	       compartment);
	padding = (align_len - (compartment_len & align_mask)) & align_mask;
	compartment_pos += compartment_len + padding;
      }
      if (sl > state->max_molecule_len) {
	state->max_molecule_len = sl;
      } else {
	if (sl < state->min_molecule_len) {
	  state->min_molecule_len = sl;
	}
      }
      /*
      matrix_text[molecules] = (char*)&raw_molecules_text[molecules_pos];
      */
      matrix_text[molecules] = molecules_pos;
      sll = (int64_t)sl + (int64_t)1;
      padding = (align_len - (sll & align_mask)) & align_mask;
      strcpy((char *)&raw_molecules_text[molecules_pos],rctnts);
      upcase(sl,(char *)&raw_molecules_text[molecules_pos],
	     (char *)&molecules_text[molecules_pos]);
      /*
      unsorted_molecules->string = (char *)&molecules_text[molecules_pos];
      */
      unsorted_molecules->string = molecules_pos;
      unsorted_molecules->m_index  = molecules;
      unsorted_molecules->c_index  = ci;
      unsorted_molecules += 1; /* Caution address arithmetic. */

      molecules_pos += (int64_t)(sll + padding);
      rctnts[sl] = ' ';

      pos += sl; /* Caution address arithmetic. */
      rctnts += sl;
      molecules += 1;
      if (side < 0) {
	reaction->num_reactants += 1;
      } else {
	reaction->num_products += 1;
      }
      skip = count_ws(rctnts);
      
      pos += skip; /* Caution address arithmetic. */
      rctnts += skip;

      if (pos < len) {
	if (rctnts[0] != '+') {
	  fprintf(stderr,"parse_reactions_file: Error, missing + "
		  "between reactants, reactions file line was \n%s\n",
		  rctnts_p);
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
  *molecules_pos_p   =  molecules_pos;
  *compartment_pos_p =  compartment_pos;
  *molecules_p       =  molecules;
  *cmpts_p           =  cmpts;
  return (success);
}
