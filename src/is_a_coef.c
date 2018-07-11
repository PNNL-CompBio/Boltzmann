/* is_a_coef.c
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
#include "is_a_coef.h"

int is_a_coef(int sl, char *line, double *coeff_p) {
  /*
    Determine whether the first sl characters of line are all digits
    As we move into the realm of floating point coefficients this routine
    becomes, determine whether the first token in line is a floating point
    coefficient.
    "abc" needs to fail
    "-3" needs to succeed
    "4" needs to succeed.
    "1.2" needs to succeed.
    "1.2e2" needs to succeed.
    "2-hdroyxy-quinone" needs to fail.
    
    Called by: count_molecules, count_molecules_and_cmpts, parse_side_line.
  */
  /*
  int i;
  int result;
  int c;
  int sp;
  result = 0;
  if (line) {
    result = 1;
    sp = 0;
    if (line[0] == '-') {
      sp = 1;
    }
    for (i=sp;i<sl;i++) {
      c = (int)line[i];
      if ((c < 48) || (c > 57)) {
	result = 0;
	break;
      } 
    }
  }
  */
  double coeff;
  int ns;
  int result;
  char c;
  char c_sl;
  c_sl = line[sl];
  line[sl] = '\0';
  /*
    Check for a - token in which case we want to return a 1 and a value of
    -1 for the coeff.
  */
  if (sl == 1) {
    if (line[0] == '-') {
      result = 1;
      *coeff_p = -1.0;
    }
  } else {
    ns = sscanf(line,"%le%c",&coeff,&c);
    result = 0;
    if (ns == 1) {
      result = 1;
      *coeff_p = coeff;
    } else {
      *coeff_p = 0.0;
    }
  }
  line[sl] = c_sl;
  return(result);
}
