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

int is_a_coef(int sl, char *line) {
  /*
    Determine whether the first sl characters of line are all digits
    Called by: count_molecules
  */
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
  return(result);
}
