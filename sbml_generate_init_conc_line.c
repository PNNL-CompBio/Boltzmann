/* sbml_generate_init_conc_line.c
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
#include "sbml_find_string.h"

#include "sbml_generate_init_conc_line.h"
void sbml_generate_init_conc_line(struct sbml2bo_struct *state,
				  int species_count,
				  char *species,
				  char *species_name,
				  char *kegg_id,
				  char *comp,
				  double multiplier,
				  double recip_volume,
				  double init_conc,
				  double init_amt,
				  int variable) {
  /*
    Print an initial concentration line for the species in the init_concs file
    and add on to the species_to_json_id map. 
    Called by: sbml_process_slist_of_species
    Calls:     boltzmannize_string, sbml_find_string,
  */
  struct t2js_struct *kg2js_data;
  struct t2js_struct *ms2js_data;
  double amt;
  double default_recip_volume;
  char vc[8];
  char bspecies_c[1024];
  char bkegg_id_c[1024];
  char bname_c[1024];
  char *specid;
  char *translation;
  char *bspecies;
  char *bname;
  char *tspecies;
  char *bkegg_id;
  char **spec_ids;
  char **translations;
  char **kegg_ids;
  char **kjs_ids;
  char **json_ids;
  char **ms_ids;
  char **mjs_ids;
  int64_t num_kegg_ids;
  int64_t num_ms_ids;

  int kegg_index;
  int not_found;

  int json_index;
  int ms_index;

  FILE *lfp;
  FILE *error_fp;
  FILE *concs_fp;
  FILE *id_name_fp;
  lfp          = state->log_fp;
  concs_fp     = state->concs_fp;
  id_name_fp   = state->id_name_fp;
  translations = state->translations;
  spec_ids     = state->spec_ids;
  kg2js_data   = state->kg2js_data;
  ms2js_data   = state->ms2js_data;
  json_ids     = state->json_ids;

  bspecies     = (char*)&bspecies_c[0];
  bkegg_id     = (char*)&bkegg_id_c[0];
  bname        = (char*)&bname_c[0];
  num_kegg_ids = kg2js_data->num_ids;
  kegg_ids     = kg2js_data->dictionary_ids;
  kjs_ids      = kg2js_data->json_ids;
  num_ms_ids   = ms2js_data->num_ids;
  ms_ids       = ms2js_data->dictionary_ids;
  mjs_ids      = ms2js_data->json_ids;

  if (lfp == NULL) {
    error_fp = stderr;
  } else {
    error_fp = lfp;
  }
  /*
    First we need to conver species,species_name,kegg_id to a json_id
    and build the mapping.
    The algorithm is roughly as follows:
     0. Copy species into the specid_2_json_strings buffer.
        pointer to the position in buffer is in spec_ids[species_count]
	Position in buffer for translation is translations[species_count]

     1. Copy species into spec_id slot 

     2. Check if we had a kegg_id show up, if so, look that up,
        and if present put the corresponding json_id into the tranlsation slot
	in specid_2_json_strings buffer. 

     3.  Boltzmannize (convert to upper case and replace spaces with underscores
         the species variable, to bspecies. 

     4.  Then check in succession for bspecies being a kegg_id, json_id
         or modelseed_id and use the corresponind json_id as a translation.

     5.  If bspecies didn't translate to a json_id,
         boltzmannize species_name to bname and check to see if bname
	 is a kegg_id, json_id or modelseed_id.

     6.  If none of the above generat a json_id, just use bname if 
         species_name existed, else don't translate, just use species.
  */
  specid = spec_ids[species_count];
  translation = translations[species_count];
  strcpy(specid,species);
  not_found = 1;
  if (kegg_id[0] != 0) {
    strcpy(bkegg_id,kegg_id);
    boltzmannize_string(bkegg_id);
    kegg_index = sbml_find_string(bkegg_id,kegg_ids,num_kegg_ids);
    if (kegg_index >= 0) {
      strcpy (translation,kjs_ids[kegg_index]);
      not_found = 0;
    } else {
      /*
	We don't have an entry for this kegg id. Log it to file.
	I want to know about these.
      */
      if (lfp) {
	strcpy (bname,species_name);
	boltzmannize_string(bname);
	fprintf(lfp,
		"kegg_id in sbml_file not found in our database: %s :: %s\n",
		bkegg_id,bname);
	fflush(lfp);
      }
    }
  }
  if (not_found) {
    /*
      Didn't have a known kegg_id attached, see if specid is a
      known json_id.
    */
    strcpy (bspecies,species);
    boltzmannize_string(bspecies);
    kegg_index = sbml_find_string(bspecies,kegg_ids,num_kegg_ids);
    if (kegg_index >= 0) {
      strcpy (translation,kjs_ids[kegg_index]);
      not_found = 0;
    }
  }
  if (not_found) {
    json_index = sbml_find_string(bspecies,json_ids,num_kegg_ids);
    if (json_index >= 0) {
      strcpy (translation,json_ids[json_index]);
      not_found = 0;
    }
  }
  if (not_found) {
    /*
      spec_id was not a json_id nor a kegg_id, check for a modelseed id.
    */
    ms_index = sbml_find_string(bspecies,ms_ids,num_ms_ids);
    if (ms_index >= 0) {
      strcpy(translation,mjs_ids[ms_index]);
      not_found = 0;
    }
  }
  if (not_found) {
    /*
      Didn't have a known kegg_id attached, and specid wasn't itself
      a kegg_id,known json_id, nor modelseed id.
      Check to see if the species_name field was set and whether or not
      it matches a kegg_id, json_id, or modelseed id.
    */
    if (species_name[0] != '\0') {
      strcpy (bname,species_name);
      boltzmannize_string(bname);
      kegg_index = sbml_find_string(bname,kegg_ids,num_kegg_ids);
      if (kegg_index >= 0) {
	strcpy (translation,kjs_ids[kegg_index]);
	not_found = 0;
      }
      if (not_found) {
	json_index = sbml_find_string(bname,json_ids,num_kegg_ids);
	if (json_index >= 0) {
	  strcpy (translation,json_ids[json_index]);
	  not_found = 0;
	}
      }
      if (not_found) {
	/*
	  species_name was not a json_id nor a kegg_id, 
	  check for a modelseed id.
	*/
	ms_index = sbml_find_string(bname,ms_ids,num_ms_ids);
	if (ms_index >= 0) {
	  strcpy(translation,mjs_ids[ms_index]);
	  not_found = 0;
	}
      }
      if (not_found) {
	/*
	  just use species name as translation
	*/
        strcpy(translation,bname);
      }
    } else {
      /*
	specid didn't match a kegg_id, json_id, nor modelseed_id,
	so just us it as itself.
      */
      strcpy(translation,bspecies);
    }
  }
  if (id_name_fp != NULL) {
    fprintf(id_name_fp,"%s\t%s",specid,translation);
    if (not_found) {
      fprintf(id_name_fp," (matching json_id not found)\n");
    } else {
      fprintf(id_name_fp,"\n");
    }
  }
  /*
    Now work on the rest of the conc_init line.
  */
  vc[0]       = 'F';
  vc[1]       = 'V';
  amt = init_amt;
  if ((amt <= 0.0) && (init_conc > 0.0)) {
    amt = init_conc;
  }
  if (amt > 0.0) {
    if (comp[0] != '\0') {
      /*
	We had a non null compartment.
      */
      amt *= multiplier * recip_volume;
      fprintf(concs_fp,"%s:%s\t%le\t%c\n",
	      translation,comp,amt,(char*)&vc[variable]);
    } else {
      /*
	Compartment did not exist.
      */
      amt *= recip_volume;
      fprintf(concs_fp,"%s\t%le\t%c\n",
	      translation,init_amt,vc[variable]);
    }
  }
}
