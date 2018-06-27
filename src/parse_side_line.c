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

#include "boltzmann_structs.h"

#include "count_ws.h"
#include "count_nws.h"
#include "upcase.h"
#include "is_a_coef.h"
#include "find_colon.h"

#include "parse_side_line.h"
int parse_side_line(char *species_p, 
		    int64_t *molecules_pos_p,
		    int64_t *compartment_pos_p,
		    int *molecules_p,
		    int *cmpts_p,
		    struct reaction_struct *reaction,
		    struct state_struct *state,
		    int side) {
  /*
    Parse a left/right line, count and record reactant molecules and 
    coefficients in the reacstion matrix.

    Called by: parse_reactions_file
    Calls:     strlen, atoi, fprint, fflush

    Variable           TMF    Description

    species_p           C*I   Pointer to the string of reactants/products.
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
  struct molecule_struct *unsorted_molecules;
  struct molecule_struct *molecule;
  struct compartment_struct *unsorted_cmpts;
  struct compartment_struct *compartment;
  struct reactions_matrix_struct *rxns_matrix;
  double *recip_coeffs;
  double dside;
  int64_t coeff;
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
  int64_t max_molecules_pos;
  int64_t max_compartment_pos;

  char    *rctnts;
  char    *compartment_text;
  char    *molecules_text;
  char    *raw_molecules_text;
  char    *compartment_name;
  char    *molecule_name;
  char    *solvent_string;

  int     token_length;
  int     colon_loc;

  int     skip;
  int     ci;

  int     padding;
  int     compartment_len;

  int     cmpts;
  int     molecules;

  int     success;
  int     molecule_name_length;

  int     max_cmpt;
  int     max_mol;

  int     cpd;
  int     mpd;

  FILE *lfp;
  FILE *efp;

  success            = 1;
  molecules_pos      = *molecules_pos_p;
  compartment_pos    = *compartment_pos_p;
  molecules          = *molecules_p;
  cmpts              = *cmpts_p;
  unsorted_molecules = state->unsorted_molecules;
  unsorted_cmpts     = state->unsorted_cmpts;
  rctnts             = species_p;
  dside              = (double)side;
  rxns_matrix        = state->reactions_matrix;
  molecules_indices  = rxns_matrix->molecules_indices;
  coefficients       = rxns_matrix->coefficients;
  recip_coeffs       = rxns_matrix->recip_coeffs;
  matrix_text        = rxns_matrix->text;

  align_len    	     = state->align_len;
  align_mask   	     = state->align_mask;
  compartment_text   = state->compartment_text;
  molecules_text     = state->molecules_text;
  raw_molecules_text = state->raw_molecules_text;
  solvent_string     = state->solvent_string;
  max_cmpt           = state->number_compartments-1;
  max_mol            = state->number_molecules-1;
  max_molecules_pos  = state->molecule_text_length;
  max_compartment_pos = state->compartment_text_length;
  pos = 0;
  len = strlen(rctnts);
  while (pos < len) {
    token_length = count_nws(rctnts);
    if (token_length > 0) {
      /*
	We have a token and thus a molecule, check to see that 
	we don't walk of the end of unsorted_molecules,
	molecules_indices, coefficients, recip_coefficients, or
	matrix_text, arrays.
      */
      if (molecules > max_mol) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"parse_side_line: Saw more molecules than "
		  "size_rxn_file counted, max_mol = %d\n",max_mol);
	  fflush(lfp);
	}
	break;
      } else {
	/* 
	   Not too many molecules.
	*/
	molecules_indices[molecules] = molecules;
	if (is_a_coef(token_length,rctnts)) {
	  /*
	    Null terminate the token.
	  */
	  rctnts[token_length] = '\0';
	  if (token_length == 1) {
	    if (rctnts[0] == '-') {
	      coeff = -1.0;
	    } else {
	      coeff = atoi(rctnts);
	    }
	  } else {
	    coeff = atoi(rctnts);
	  }
	  if (side < 0) {
	    coeff = - coeff;
	  } 
	  coefficients[molecules] = coeff;
	  recip_coeffs[molecules] = 1.0/((double)coeff);
	  /*
	    Replace the null with  a space.
	  */
	  rctnts[token_length] = ' ';

	  pos += (int64_t)token_length;
	  rctnts += token_length; /* Caution address arithmetic. */
	  /*
	    Skip over the next block of white space.
	  */
	  skip   =  count_ws(rctnts);
	  
	  pos += skip;
	  rctnts += skip; /* Caution address arithmetic. */
	  
	  token_length = count_nws(rctnts);
	} else {
	  coefficients[molecules] = (int64_t)side;
	  recip_coeffs[molecules] = (double)side;
	}
	/*
	  At this juncture we need to check for
	  a "local" compartment specification of the form 
	  :compartment trailing the molecule name. We do not
	  allow spaces on either side of the semicolon.
	*/
	ci = -1;
	/*
	  Null terminate the token.
	*/
	rctnts[token_length] = '\0';
	molecule_name_length = token_length;
	colon_loc = find_colon(rctnts);
	compartment_len = 0;
	if (colon_loc >= 0) {
	  /* 
	     We had a local :compartment attached.
	     determine its length, store it and shorten the
	     length for the molecule - do not forget to 
	     allow for the terminating null.
	  */
	  compartment_len = token_length - colon_loc - 1;
	  molecule_name_length = colon_loc;
	  /*
	    Null terminate the molecle name length replacing the
	    : with a null.
	  */
	  rctnts[molecule_name_length] = '\0';
	  /*
	    check for longest and shorted compartment_name.
	  */
	  if (compartment_len > state->max_compartment_len) {
	    state->max_compartment_len = compartment_len;
	  } else {
	    if (compartment_len < state->min_compartment_len) {
	      state->min_compartment_len = compartment_len;
	    }
	  }
	  /*
	    Check to see that we will not walk of end of 
	    unworted_cmpts array.
	  */
	  if (cmpts > max_cmpt) {
	    success = 0;
	    if (lfp) {
	      fprintf(lfp,"parse_size_Line: error parsed more compartments"
		      " than counted in size_rxns_file, max_cmpt = %d\n",
		      max_cmpt);
	      fflush(lfp);
	    }
	    break;
	  } else {
	    /*
	      Store the comparment information.
	    */
	    compartment = (struct compartment_struct *)&unsorted_cmpts[cmpts];
	    compartment->string = compartment_pos;
	    compartment->volume = 0.0;
	    compartment_name = (char *)&compartment_text[compartment_pos];
	    compartment->c_index  = cmpts;
	    /*
	      Allow for the terminating \0
	    */
	    padding = (align_len - ((compartment_len+1) & align_mask)) & align_mask;
	    cpd = compartment_len + 1 + padding;
	    ci = cmpts;
	    cmpts += 1;
	    /*
	      Realy before this strcpy we should check that
	      compartment_pos + compartment_len + 1 + padding does not
	      execeed compartment_text_length allocated in 
	      alloc2_a.c
	    */
	    if ((compartment_pos + cpd) > max_compartment_pos) {
	      success = 0;
	      if (lfp) {
		fprintf(lfp,"parse_size_line: Error compartment names plus"
			" alignment padding longer than %ld measured in "
			"size_rxns_file.c\n",max_compartment_pos);
		fflush(lfp);
			
	      }
	      break;
	    }
	    strcpy(compartment_name,(char*)&rctnts[colon_loc+1]);
	    upcase(compartment_len,compartment_name);
	    compartment_pos += (int64_t)cpd;
	  }
	} /* end if colon_loc > 0 */
	/*
	  Check for longest and shortest molecule_name.
	*/
	if (molecule_name_length > state->max_molecule_len) {
	  state->max_molecule_len = molecule_name_length;
	} else {
	  if (molecule_name_length < state->min_molecule_len) {
	    state->min_molecule_len = molecule_name_length;
	  }
	}
	/*
	 matrix_text[molecules] = (char*)&raw_molecules_text[molecules_pos];
	*/
	matrix_text[molecules] = molecules_pos;
	sll = (int64_t)molecule_name_length + (int64_t)1;
	padding = (align_len - (sll & align_mask)) & align_mask;
	molecule_name = (char *)&molecules_text[molecules_pos];
	/*
	  Really before this strcpy we should check that
	  molecules_pos + molecule_name_length + 1 + padding,
	  does not exceed the allocated molecule_text_length
	  in alloc2_a.c
	*/
	mpd = sll + padding;
	if ((molecules_pos + mpd) > max_molecules_pos) {
	  success = 0;
	  if (lfp) {
	    fprintf(lfp,"parse_sideline: Error molecule name lengths plus"
		    " alignment padding longer than measured in "
		    "size_rxns_file.c, = %ld\n",max_molecules_pos);
	    fflush(lfp);
	  }
	} else {
	  strcpy(molecule_name,rctnts);
	  upcase(molecule_name_length,molecule_name);
	}
	/*
	  unsorted_molecules->string = (char *)&molecules_text[molecules_pos];
	*/
	molecule = (struct molecule_struct *)&unsorted_molecules[molecules];
	molecule->string = molecules_pos;
	molecule->m_index  = molecules;
	molecule->c_index  = ci;
	if (strcmp(molecule_name,solvent_string) == 0) {
	  molecule->solvent = 1;
	} else {
	  molecule->solvent = 0;
	}
	/*
	  Set the variable field to -1 to track prescence of molecule
	  in the intial concentrations file.
	*/
	molecule->variable = -1;

	molecules_pos += (int64_t)mpd;
	if (compartment_len > 0) {
	  rctnts[molecule_name_length] = ':';
	}
	if (side < 0) {
	  reaction->num_reactants += 1;
	} else {
	  reaction->num_products += 1;
	}
      } /* end else not too many molecules.n */
      rctnts[token_length] = ' ';

      pos += token_length; 
      rctnts += token_length; /* Caution address arithmetic. */
      molecules += 1; 
      skip = count_ws(rctnts);
      
      pos += skip; 
      rctnts += skip; /* Caution address arithmetic. */

      /*
	If there are more tokens on the line check to see that
	the next tokein is + or a - and skip over it.
      */
      if (pos < len) {
	/*
	  If connecting sign is a + skip over it.
	*/
	if (rctnts[0] == '+')  {
	  pos += 1;
	  rctnts += 1; /* Caution address arithmetic.*/

	  skip = count_ws(rctnts);
	
	  pos  += skip;
	  rctnts += skip; /* Caution address arithmetic.*/
	} else {
	  /*
	    If connecting sign is a minus, don't skip over it to 
	    allow it to be parsed as a minus 1 stoichiometric coefficient
	    next pass through the loop.
	  */
	  if (rctnts[0] != '-') {
	    success = 0;
	    if (lfp) {
	      fprintf(lfp,"parse_side_line: Error, missing + or - "
		      "between species, reactions file line was \n%s\n",
		      species_p);
	      fflush(lfp);
	    }
	    break;
	  } 
	}
      }
    } /* end if (token_length > 0)  - a reactant was found . */
  } /* end while (pos > len) */
  *molecules_pos_p   =  molecules_pos;
  *compartment_pos_p =  compartment_pos;
  *molecules_p       =  molecules;
  *cmpts_p           =  cmpts;
  return (success);
}
