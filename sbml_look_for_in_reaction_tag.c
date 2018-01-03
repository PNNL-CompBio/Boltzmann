/* sbml_look_for_in_reaction_tag.c
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

#include "sbml_look_for_in_reaction_tag.h"

int sbml_look_for_in_reaction_tag (char *line) {
/*
  Identify tags found within the reaction tag, including itself.
  Called by: sbml_process_list_of_reactions.
  Calls:     strncmp
tag val      tag,length

below is a list of the tags in sorted alpha ordering.

returned
 32  	   "</listOfProducts>",17
  8  	   "</listOfReactants>",18
256  	   "</reaction>",11

 16  	   "<listOfProducts",15
  4  	   "<listOfReactants",16

128  	   "<kineticLaw>",12
  2  	   "<reaction",9
 64  	   "<speciesReference",17
  1        none of the above.
*/
  int tag;
  int padi;
  tag = 1;
  if (strlen(line) > 8) {
    if (line[0] == '<') {
      /*
	We go the first character correct.
      */
      if (line[1] == '/') {
	if (line[2] == 'l') {
	  if (strncmp(line,"</listOfProducts>",17) == 0) {
	    tag = 32;
	  } else {
	    if (strncmp(line,"</listOfReactants>",18) == 0) {
	      tag = 8;
	    }
	  }
	} else {
	  if (line[2] == 'r') {
	    if (strncmp(line,"</reaction>",11) == 0) {
	      tag = 256;
	    }
	  }
	}
      } else {
	if (line[1] == 'l') {
	  if (strncmp(line,"<listOfProducts",15) == 0) {
	    tag = 16;
	  } else {
	    if (strncmp(line,"<listOfReactants",16) == 0) {
	      tag = 4;
	    }
	  }
	} else {
	  if (line[1] == 'k') {
	    if (strncmp(line,"<kineticLaw>",12) == 0) {
	      tag = 128;
	    }
	  } else {
	    if (line[1] == 'r') {
	      if (strncmp(line,"<reaction",9) == 0) {
		tag = 2;
	      }
	    } else {
	      if (line[1] == 's') {
		if (strncmp(line,"<speciesReference",17) == 0) {
		  tag = 64;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  return(tag);
}
