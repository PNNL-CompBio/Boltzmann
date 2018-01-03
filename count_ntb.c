/* count_ntb.c
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

#include "count_ntb.h"

int count_ntb(char *line) {
  /*
    Return the number of leading non-tab characters ('\09') 
    in the character string, line.
    Called by: read_kg2js, read_ms2js
    Calls      strlen(intrinsic).
  */
  int ntb_chars;
  int i;
  int c;
  int tab;
  /*
    Want to return 0 if line has zero length.
  */
  tab = 9;
  ntb_chars = 0;
  for (i=0;i<strlen(line);i++) {
    c  = (int)line[i];
    if (c == tab) break;
    ntb_chars += 1;
  }
  return (ntb_chars);
}
