/* count_molecules.c
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

#include "count_ws.h"
#include "count_nws.h"
#include "is_a_coef.h"
#include "count_molecules.h"

int count_molecules(char* molecules_line, int64_t *molecules_len) {

  /*
    Count the number of molecules occuring in a reactants (LEFT) or
    products (RIGHT) line, assumes LEFT or RIGHT prefix has already
    been remove from molecules_line.
    Called by size_rxns_file
  */
  char *line;
  int line_len;
  int pos;

  int cur_pos;
  int c;

  int molecules;
  int sl;
  int wsl;

  line_len = strlen(molecules_line);
  molecules = 0;
  line     = molecules_line;
  if (line_len > 0) {
    /*
      Skip over leading white space.
    */
    wsl = count_ws(line);
    pos = wsl;
    line += wsl;
    while (pos < line_len) {
      /*
	Get length of molecules name.
      */
      sl = count_nws(line);
      if (sl > 0) {
	molecules += 1;
	pos += sl;
	if (is_a_coef(sl,line)) {
	  /* now it may be that a molecule is preceded by a
	     coefficient
	     if that is the case then we need to skip the coefficient
	     and the following white space.
	  */
	  line += sl;
	  wsl  = count_ws(line);
	  pos += wsl;
	  line += wsl;
	  sl   = count_nws(line);
	  pos += sl;
	}
	line += sl;
      }
      *molecules_len += (int64_t)sl;
      /*
	chew up white space after molecule.
      */
      wsl   = count_ws(line);
      pos += wsl;
      line  += wsl;
      if (pos < line_len) {
	/*
	  if we are not at the end of the line, 
	  the next character needs to be a '+'
	*/
	if (line[0] != '+') {
	  fprintf(stderr,"count_molecules: Error - malformed line\n"
		  "expecting a +, got a %c\n",line[0]);
	  fflush(stderr);
	}
	pos += 1;
	line += 1;
      }
      wsl = count_ws(line);
      pos += wsl;
      line += wsl;
    } /* end while (pos ...) */
  } /* end if (line ...) */
  return(molecules);
}
