/* sbml_process_list_of_products_tag.c
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
#include "sbml_read_key_value.h"

#include "sbml_process_list_of_products_tag.h"

void sbml_process_list_of_products_tag(int max_key_len,
				       int max_val_len,
				       char *key,
				       char *value,
				       char *buffer,
				       char *comp,
				       int  not_in_tag,
				       int  *enclosing_tag_p,
				       int  *species_count_p,
				       FILE *rxns_fp,
				       FILE *error_fp) {
  /*
    Process the <listOfProducts tag setting up a right compartment
    if it is specified and printent a "RIGHT" line in the rxns.dat file.
    Called by: sbml_process_list_of_reactions.
    Calls:     count_ws,count_nws,sbml_read_key_value
               strlen,strcmp,strcpy,fprintf,fflush
  */
  char *line;
  int enclosing_tag;
  int species_count;
  int nb;
  int ns;
  int tl;
  int padi;
  line          = buffer;
  species_count = *species_count_p;
  enclosing_tag = *enclosing_tag_p;
  /*
    Skip over the <listOfProducts tag.
  */				   
  line += 15;
  nb = count_ws(line);
  line += nb; /* Caution address arithmetic */
  ns = strlen(line);
  while (ns != 0)  {
    if (line[0] == '>') {
      /*
	End of list_of_products tab,
	print a RIGHT_COMPARTMENT line if compartment is not null,
	Then print the start of the RIGHT line.
      */
      if (comp[0] != '\0') {
	fprintf(rxns_fp,"RIGHT_COMPARTMENT\t%s\n",comp);
	comp[0] = '\0';
      }
      fprintf(rxns_fp,"RIGHT\t");
      enclosing_tag = not_in_tag;
      species_count = 0;
      ns = 0;
    } else {
      tl = sbml_read_key_value(line,key,value,max_key_len,max_val_len);
      if (tl <= 0) {
	ns = 0;
      } else {
	line += tl;
	nb = count_ws(line);
	line += nb; /* Caution address arithmetic */
	ns = strlen(line);
	if (strcmp(key,"compartment") == 0) {
	  strcpy(comp,value);
	  nb = count_ws(line);
	  line += nb; /* Caution address arithmetic */
	  ns = count_nws(line);
	} else {
	  /*
	  fprintf(error_fp,"sbml_process_list_of_reactions: Error "
		  "Unrecognized key=value pair in listOfProducts "
		  "tag:%s=%s\n",key,value);
	  fflush(error_fp);
	  */
	}
      } /* end else found a key=value field in the listOfReactans tag.*/
    } /* end else didn't find end of listOfProducts tag */	  
  } /* end while (ns != 0) */
  *enclosing_tag_p = enclosing_tag;
  *species_count_p = species_count;
}
