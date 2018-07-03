/* sbml_process_species_reference_tag.c
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
#include "sbml_lookup_speciesref_attribute.h"
#include "sbml_find_string.h"

#include "sbml_process_species_reference_tag.h"
void sbml_process_species_reference_tag(struct sbml2bo_struct *state,
					int max_key_len,
					int max_val_len,
					char *key,
					char *value,
					char *buffer,
					char *species,
					char *comp, 
					int not_in_tag, 
					int *species_count_p,
					double *coefficient_p,
					int *enclosing_tag_p,
					FILE *rxns_fp,
					FILE *error_fp) {

  /*
    Use the json_id translation of the species tag.
    Called by: sbml_process_list_of_reactions.
    Calls:     count_ws,count_nws,sbml_read_key_value,
               sbml_lookup_speciesref_attribute,
	       sbml_find_string,
               strncmp,sscanf,fprintf,fflush
  */
  char **spec_ids;
  char **translations;
  char *line;
  double coefficient;
  int nb;
  int ns;

  int species_count;
  int num_species;

  int enclosing_tag;
  int tl;

  int tag;
  int trans_index;


  num_species   = state->num_species;
  spec_ids     = state->spec_ids;
  translations = state->translations;
  line          = buffer;
  species_count = *species_count_p;
  coefficient   = *coefficient_p;
  enclosing_tag = *enclosing_tag_p;
  line += 17;
  nb = count_ws(line);
  line += nb; /* Caution address arithmetic */
  ns = count_nws(line);
  while (ns != 0) {
    if ((strncmp(line,"/>",2) == 0) || (line[0] == '>')){
      /*
	We have reached the end of the speciesReference tag,
	print out the stoichiometry, species name and compartment
	if it exists, and preceding + sign if species_count is > 0.
      */
      if (species_count > 0) {
	fprintf(rxns_fp,"\t+\t");
      }
      if (coefficient != 1.0) {
	fprintf(rxns_fp,"%le ",coefficient);
      }
      fprintf(rxns_fp,"%s",species);
      if (comp[0] != '\0') {
	fprintf(rxns_fp,":%s",comp);
	comp[0] = '\0';
      }
      /*
      coefficient = 1.0;
      */
      species_count += 1;
      enclosing_tag = not_in_tag;
      ns = 0;
    } else {
      /*
	look for key=value triples, expecting keys of
	"species", "stoichiometry" "constant" and possibly "compartment"
      */
      tl = sbml_read_key_value(line,key,value,max_key_len,max_val_len);
      if (tl <= 0) {
	ns = 0;
      } else {
	line += tl;
	nb = count_ws(line);
	line += nb; /* Caution address arithmetic */
	tag = sbml_lookup_speciesref_attribute(key);
	
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
	      constant, ignored for now.
	    */
	    break;
	  case 2:
	    /*
	      species
	    */
	    strcpy(species,value);
	    trans_index = sbml_find_string(species,spec_ids,num_species);
	    if (trans_index >= 0) {
	      strcpy(species,translations[trans_index]);
	    }
	    break;
	  case 3:
	    /*
	      stoichiometry.
	    */
	    sscanf(value,"%le",&coefficient);
	    break;
	  } /* end switch */  
	} else {
	  fprintf(error_fp,"sbml_process_species_reference_tag: "
		  "Error found unexpected speciesref key :%s\n",key);
	  fflush(error_fp);
	} 
      } /* end else we had key=value triple */
    } /* end else not at end of speciesReference tag. */
  } /* end while (ns ...) */
  *species_count_p  = species_count;
  *coefficient_p    = coefficient;
  *enclosing_tag_p  = enclosing_tag;
}
