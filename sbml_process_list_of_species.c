/* sbml_process_list_of_species.c
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
#include "sbml_start_species_def.h"
#include "sbml_parse_species_key_value.h"
#include "sbml_generate_init_conc_line.h"

#include "sbml_process_list_of_species.h"
int sbml_process_list_of_species(FILE *sbml_fp,
				 char *sbml_buffer,
				 int sbml_buffer_len, 
				 struct sbml2bo_struct *state) {
  /*
    Process the <listOfSpecies> tag contents, and generate
    the concs.in file.
    Called by: parse_sbml
    Calls:     fprintf, feof, fgets, fflush, strcpy, strcmp, strncmp,
               count_ws, count_nws,
	       sbml_start_species_def,
	       sbml_parse_species_key_value,
               sbml_generate_init_conc_line
  */
  struct compartment_struct *compartment;
  struct t2js_struct *ms2js_data;
  char comp_c[1024];
  char species_c[1024];
  char species_name_c[1024];
  char kegg_id_c[1024];
  char *kegg_id;
  char *comp;
  char *species;
  char *species_name;
  char bcf[8];
  char *concs_file;
  char *bc;
  char *line;
  
  double init_amt;
  double init_conc;
  double multiplier;
  double recip_volume;
  double default_recip_volume;

  int variable;
  int in_species_tag;

  int in_list_of_species;
  int success;

  int nb;
  int ns;

  int in_species_group;
  int i;

  int tl;
  int species_count;
  
  FILE *lfp;
  FILE *error_fp;
  FILE *concs_fp;
  FILE *id_name_fp;

  ms2js_data  = state->ms2js_data;
  success     = 1;
  comp        = (char *)&(comp_c[0]);
  species     = (char *)&(species_c[0]);
  species_name = (char *)&(species_name_c[0]);
  kegg_id     = (char *)&(kegg_id_c[0]);
  bc          = (char *)&bcf[0];

  lfp         = state->log_fp;
  concs_fp    = state->concs_fp;
  id_name_fp  = state->id_name_fp;
  if (lfp == NULL) {
    error_fp = stderr;
  } else {
    error_fp = lfp;
  }
  /*
    For now use the default volume for the volume line, we will 
    read in the individual compartment volumes in the
    listOfCompartments section.
  */
  default_recip_volume = 1.0e15;
  fprintf(concs_fp,"VOLUME 1.0e-15\n");
  /*
    At some point in time we will want to compute the units field
    from the substanceUnits field but for now us 1.e-9
  */
  fprintf(concs_fp,"CONC_UNITS 1.0e-9\n");
  in_list_of_species = 1;
  in_species_tag = 0;
  in_species_group = 0;
  species_count = 0;
  while (!feof(sbml_fp) && (in_list_of_species == 1)) {
    /* 
       Now we are in the listOfSpecies section.    
    */
    fgets(sbml_buffer,sbml_buffer_len,sbml_fp);
    line = sbml_buffer;
    nb = count_ws(line);
    line += nb; /* Caution address arithmetic */
    /*
      First check for an </listOfSpecies> tag.
    */
    if (strncmp(line,"</listOfSpecies>",16) == 0) {
      in_list_of_species = 0;
    } else {
      if (in_species_group == 0) {
	if (in_species_tag == 0) {
	  if (strncmp(line,"<species",8) == 0) {
	    /*
	      Initialize species attributes and set in_species_tag
	      and in_species group to 1.
	    */
	    sbml_start_species_def(species,species_name,
				   kegg_id,comp,
				   default_recip_volume,
				   &multiplier,
				   &recip_volume,
				   &init_conc,
				   &init_amt,
				   &variable,
				   &in_species_tag,
				   &in_species_group);
	    line +=8; /* Caution address arithmetic */
	    /*
	      Skip over white space.
	    */
	    nb = count_ws(line);
	    line += nb; /* Caution address arithmetic */
	  } /* end if <species line */
	} /* end if not in <species tag */
      } /* end if not in <species></species> group */
      if (in_species_tag) {
        ns = count_nws(line);
        while (ns != 0) {
          /*
            check for end tag.
          */
          if ((strncmp(line,"/>",2)==0) || (line[0] == '>')) {
            ns = 0;
            in_species_tag = 0;
          } else {
            /*
	      Not the end of species  tag.
            */
	    tl = sbml_parse_species_key_value(state,line,
					      default_recip_volume,
					      bc,comp,
					      species,species_name,
					      &recip_volume,&multiplier,
					      &init_amt,&init_conc,
					      &variable);
            if (tl <= 0) {
	      /*
		we have reached the end of this line in the file,
		but not yet the closing tag of <species .
	      */
              ns = 0;
            } else {
              line += tl;
            } /* end else we found a keyword=value triple. */
          } /* end else we did not find an end of the species tag */
        } /* end while (ns != 0) */
      } else {
	/*
	  Not in species tag.
	  Check for an end of species group tag, "</species>".
	*/
	if (strncmp(line,"</species>",10) == 0) {
	  /*
	    Generate species line and translation map.
	  */
	  sbml_generate_init_conc_line(state,species_count,
				       species,species_name,
				       kegg_id,comp,
				       multiplier,
				       recip_volume,
				       init_conc,init_amt,
				       variable);
	  species_count += 1;
	  in_species_group = 0;
	} else {
	  /*
	    Not at end of <species></species> group
	    look for a kegg_id if we have not gotten one.;
	  */
	  if (kegg_id[0] == '\0') {
	    if (strncmp(line,"<rdf:li rdf:resource=\"http://identifiers.org/kegg.compound/",59)==0) {
	      line += 59;
	      sscanf(line,"%s",kegg_id);
	      for (i=0;i<strlen(kegg_id);i++) {
		if (kegg_id[i] == '"') {
		  kegg_id[i] = '\0';
		  break;
		}
	      }
	    }
	  } /* end if no kegg_id yet */
	} /* end else not at end of <species></species> group */
      } /* end else not in species tag. */
    } /* end else not at </listOfSpecies> tag */ 
  } /* end while (!feof...) */
  return (success);
}
