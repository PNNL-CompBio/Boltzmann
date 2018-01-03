/* boltzmannize_string.c
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

#include "boltzmannize_string.h"
void boltzmannize_string(char *string) {
  int i;
  int len;
  int ic;
  int j;
  char *spec_char;
  len = strlen(string);
  j = 0;
  for (i=0;i<len;i++) {
    ic = (int)string[i];
    /*
      Use any form of white space as a string terminator.
    */
    if (ic <= 32) {
      string[j] = '_';
      j += 1;
    } else {
      if ((ic > 96) && (ic < 123)) {
	string[j] = ic - 32;
	j += 1;
      } else {
	/*
	  Remove double quotes.
	*/
	if (ic != 34) {
	  string[j] = ic;
	  j += 1;
	}
      }
    }
  } /* end for ... */
  string[j] = '\0';
  /*
    Now that it has been upper cased  we need to check for the 
    &*; special html characters. For now the only two that we have seen
    are &GT; and &APOS
  */
  j=0;
  len = strlen(string);
  for (i=0;i<len;) {
    ic = (int)string[i];
    if (ic == 38) {
      spec_char = (char *)&string[i];
      if (strncmp(spec_char,"&GT;",4) == 0) {
	string[j] = '>';
	i += 4;
	j += 1;
      } else {
	if (strncmp(spec_char,"&APOS;",6) == 0) {
	  string[j] = '\'';
	  i += 6;
	  j += 1;
	} else {
	  string[j] = string[i];
	  i += 1;
	  j += 1;
	}
      }
    } else {
      string[j] = string[i];
      i+= 1;
      j+= 1;
    }
  }
  string[j] = '\0';
}
