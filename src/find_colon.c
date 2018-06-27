/* find_colon.c
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
#include "find_colon.h"
int find_colon(char *rctnts) {
  /*
    Called by: parse_side_line, translate_regulation_metabolites
    Return the position of a the leftmost colon in a string.
    Returns -1 if no colon was found.
  */
  int i;
  int colon_loc;
  colon_loc = -1;
  for (i=0;i<strlen(rctnts);i++) {
    if (rctnts[i] == ':') {
      colon_loc = i;
      break;
    }
  }
  return(colon_loc);
}
