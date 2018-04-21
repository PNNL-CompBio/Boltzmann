/* sbml_process_reaction_tag.c
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

#include "sbml_read_key_value.h"
#include "sbml_lookup_reaction_attribute.h"
#include "count_ws.h"
#include "count_nws.h"

#include "sbml_process_reaction_tag.h"

void sbml_process_reaction_tag(int max_key_len,
			       int max_val_len,
			       char *key,
			       char *value,
			       char *buffer,
			       char *comp,
			       char *units,
			       char *rxn_id,
			       char *rxn_name,
			       double *delta_g0_p,
			       FILE *rxns_fp,
			       FILE *error_fp) {
  /*
    Process the fields of an sbml <reaction> tag.
    Called by: sbml_process_list_of_reactions
    Calls:     sbml_read_key_value,
               sbml_lookup_reaction_attribute,
	       count_ws,
	       count_nws,
	       strcpy, sscanf, fprintf, fflush
	       
  */
  double delta_g0;
  char   *line;
  int nb;
  int ns;
  int tl;
  int tag;
  /*
    "<reaction" tag has already been skipped over. 
    Look for end of reaction tag,">"
  */
  line = buffer;
  nb = count_ws(line);
  line += nb; /* Caution address arithmetic */
  ns = count_nws(line);
  delta_g0  = *delta_g0_p;
  rxn_name[0] = '\0';
  while (ns != 0) {
    if (line[0] == '>') {
      /*
	We have hit the end of the reaction tag,
	generate  the reaction label, and
	compartment if there was one, reset
	the compartment, and break out of this loop.
      */
      if (rxn_name[0] == '\0') {
	fprintf(rxns_fp,"REACTION\t%s\n",rxn_id);
      } else {
	fprintf(rxns_fp,"REACTION\t%s :: %s\n",rxn_id,rxn_name);
      }
      if (comp[0] != '\0') {
	fprintf(rxns_fp,"COMPARTEMNT\t%s\n",comp);
      }
      comp[0] = '\0';
      ns = 0;
    } else {
      /*
	look for key=value triples, expecting keys of
	"id", "reversible" "fast" and possibly "compartment"
      */
      tl = sbml_read_key_value(line,key,value,max_key_len,max_val_len);
      if (tl <= 0) {
	ns = 0;
      } else {
	line += tl;
	tag = sbml_lookup_reaction_attribute(key);
	nb = count_ws(line);
	line += nb; /* Caution address arithmetic */
	ns = count_nws(line);
	if (tag >= 0) {
	  switch (tag) {
	  case 0:
	    /*
	      compartment.
	    */
	    strcpy (comp,value);
	    break;
	  case 1:
	    /*
	      delta_g0
	    */
	    sscanf(value,"%le",&delta_g0);
	    break;
	  case 2:
	    /*
	      delta_g0_units
	    */
	    strcpy(units,value);
	  case 3:
	    /*
	      fast, ignored for now.
	    */
	    break;
	  case 4:
	    /*
	      id
	    */
	    strcpy(rxn_id,value);
	    break;
	  case 5:
	    /*
	      reversible, ignored for now.
	    */
	    break;
	  case 6:
	    strcpy(rxn_name,value);
	    break;
	  } /* end switch */  
	} else {
	  /* Do not want to know about this for now, just ignore it.
	  fprintf(error_fp,"sbml_process_reaction_tag: "
		  "Error found unexpected reaction key :%s\n",key);
	  fflush(error_fp);
	  */
	}
      } /* end else we found a keyword=value triple. */
    } /* end else we did not find and en reaction_tag tag */
  } /* end while (ns != 0), in reaction tag. */
  *delta_g0_p = delta_g0;
}
