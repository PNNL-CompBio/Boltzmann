/* count_molecules_and_cmpts.c
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
#include "is_a_coef.h"
#include "find_colon.h"
#include "count_molecules_and_cmpts.h"

void count_molecules_and_cmpts(char* molecules_line, int *num_molecules_p, int *num_compartments_p, int64_t *molecules_len_p, int64_t *compartment_len_p, FILE *lfp) {

  /*
    Count the number of molecules and compartments 
    occuring in a reactants (LEFT) or
    products (RIGHT) line, assumes LEFT or RIGHT prefix has already
    been removed from molecules_line. Also count their lengths.
    the contents of the num_molecules_p, num_compartments_p, and
    molecules_len_p and compartment_lenp are incremented by the counts and
    lengths computed here.
    Called by: size_rxns_file
  */
  double coef;
  char *line;
  int64_t molecules_len;
  int64_t compartment_len;
  int line_len;
  int pos;

  int molecules;
  int compartments;

  int token_length;
  int molecule_name_length;

  int wsl;
  int colon_loc;

  line_len = strlen(molecules_line);
  molecules = 0;
  compartments = 0;
  molecules_len = (int64_t)0;
  compartment_len = (int64_t)0;
  line     = molecules_line;
  if (line_len > 0) {
    /*
      Skip over leading white space.
    */
    wsl = count_ws(line);
    pos = wsl;
    line += wsl; /* Caution address arithmetic here. */
    while (pos < line_len) {
      /*
	Get length of molecules name.
      */
      token_length = count_nws(line);
      if (token_length > 0) {
	molecules += 1;
	pos += token_length;
	if (is_a_coef(token_length,line,&coef)) {
	  /* now it may be that a molecule is preceded by a
	     coefficient
	     if that is the case then we need to skip the coefficient
	     and the following white space, and get the next token.
	  */
	  line += token_length; /* Caution address arithmetic here. */
	  wsl  = count_ws(line);
	  pos += wsl;
	  line += wsl; /* Caution address arithmetic here. */
	  token_length   = count_nws(line);
	  pos += token_length;
	}
	/*
	  Set the first white_space after the token to null so it becomes a 
	  string to pass to find_colon.
	*/
	line[token_length] = '\0';
	colon_loc = find_colon(line);
	if (colon_loc >= 0) {
	  compartments += 1;
	  molecule_name_length = colon_loc;
	  compartment_len += (token_length - colon_loc - 1);
	} else {
	  molecule_name_length = token_length;
	}
	line[token_length] = ' ';
	molecules_len += molecule_name_length;
      }
      line += token_length; /* Caution address arithmetic here. */
      /*
	chew up white space after molecule.
      */
      wsl   = count_ws(line);
      pos += wsl;
      line  += wsl; /* Caution address arithmetic here. */
      if (pos < line_len) {
	/*
	  if we are not at the end of the line, 
	  the next character needs to be a '+'
	*/
	if ((line[0] != '+') && (line[0] != '-')) {
	  if (lfp) {
	    fprintf(lfp,"count_molecules_and_cmpts: Error - malformed line\n"
		    "expecting a + or -, got a %c\n",line[0]);
	    fflush(lfp);
	  }
	}
	pos += 1;
	line += 1; /* Caution address arithmetic here. */
      }
      wsl = count_ws(line);
      pos += wsl;
      line += wsl; /* Caution address arithmetic here. */
    } /* end while (pos ...) */
    
  } /* end if (linelen > 0 ...) */
  *num_molecules_p += molecules;
  *num_compartments_p += compartments;
  *molecules_len_p += molecules_len;
  *compartment_len_p += compartment_len;
}
