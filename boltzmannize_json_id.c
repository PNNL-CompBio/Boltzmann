/* boltzmannize_json_id.c
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

#include "boltzmannize_json_id.h"
int boltzmannize_json_id(char *json_id) {
  int i;
  int len;
  int ic;
  int j;
  len = strlen(json_id);
  j = 0;
  for (i=0;i<len;i++) {
    ic = (int)json_id[i];
    if (ic <= 32) {
      json_id[j] = '_';
      j += 1;
    } else {
      if ((ic > 96) && (ic < 123)) {
	json_id[j] = ic - 32;
	j += 1;
      } else {
	/*
	  Remove double quotes.
	*/
	if (ic != 34) {
	  json_id[j] = ic;
	  j += 1;
	}
      }
    }
  }
  json_id[j] = '\0';
}
