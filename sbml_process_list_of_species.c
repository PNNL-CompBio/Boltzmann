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
#include "sbml_key_value.h"
#include "sbml_lookup_species_attribute.h"
#include "compartment_lookup.h"
#include "translate_ms_2_js.h"

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
	       sbml_key_value,
               sbml_lookup_species_attribute
  */
  struct state_struct pseudo_state;
  struct compartment_struct *compartment;
  struct ms2js_struct *ms2js_data;
  char comp_c[1024];
  char species_c[1024];
  char units_c[1024];
  char *key_c[1024];
  char *value_c[1024];
  char *species_name_c[1024];
  char *key;
  char *value;
  char *comp;
  char *species;
  char *tspecies;
  char *species_name;
  char *units;
  char bcf[8];
  char vc[8];
  char *concs_file;
  char *bc;
  char *line;
  
  double init_amt;
  double mole;
  double millimole;
  double micromole;
  double nanomole;
  double picomole;
  double femtomole;
  double multiplier;
  double recip_volume;
  double default_recip_volume;

  int variable;
  int in_species_tag;

  int in_list_of_species;
  int is_a_conc;
  
  int success;
  int tl;

  int max_key_len;
  int max_val_len;

  int nb;
  int ns;

  int tag;
  int ci;

  int translate_from_modelseed;
  int padi;
   

  FILE *lfp;
  FILE *error_fp;
  FILE *concs_fp;
  FILE *extra_fp;
  FILE *id_name_fp;

  pseudo_state.compartment_text = state->compartment_text;
  pseudo_state.sorted_cmpts = state->sorted_compartments;
  pseudo_state.nunique_compartments = state->num_cmpts;
  ms2js_data  = state->ms2js_data;
  success     = 1;
  key         = (char *)&(key_c[0]);
  value       = (char *)&(value_c[0]);
  comp        = (char *)&(comp_c[0]);
  species     = (char *)&(species_c[0]);
  species_name = (char *)&(species_name_c[0]);
  units       = (char *)&(units_c[0]);
  bc          = (char *)&bcf[0];
  vc[0]       = 'F';
  vc[1]       = 'V';
  max_key_len = 1023;
  max_val_len = 1023;
  mole = 1.0;
  millimole = 1.e-3;
  micromole = 1.e-6;
  nanomole  = 1.e-9;
  picomole  = 1.e-12;
  femtomole = 1.e-15;
  /*
    For now we hardware translate_from_modelseed to be 1,
    assuming the sbml files are using the modelseed id's.
    Could become a parameter later.
  */
  translate_from_modelseed = 1;

  lfp         = state->log_fp;
  concs_fp    = state->concs_fp;
  id_name_fp  = state->id_name_fp;
  if (lfp == NULL) {
    error_fp = stderr;
  } else {
    error_fp = lfp;
  }
  if (success) {
    /*
      For now use the default volume for the volume line, we will 
      read in the individual compartment volumes in the
      listOfCompartments section.
    */
    fprintf(concs_fp,"VOLUME 1.0e-15\n");
    /*
      At some point in time we will want to compute the units field
      from the substanceUnits field but for now us 1.e-9
    */
    fprintf(concs_fp,"CONC_UNITS 1.0e-9\n");
    comp[0] = '\0';
    species[0] = '\0';
    species_name[0] = '\0';
    init_amt=0;
    units[0] = '\0';
    variable     = 1;
    in_species_tag = 0;
    in_list_of_species = 1;
    while (!feof(sbml_fp) && (in_list_of_species == 1)) {
      /* 
	 Now we are in the listOfSpecies section.    
      */
      fgets(sbml_buffer,sbml_buffer_len,sbml_fp);
      line = sbml_buffer;
      if (in_species_tag == 0) {
	nb = count_ws(line);
	line += nb; /* Caution address arithmetic */
	/*
	  First check for an </listOfSpecies> tag.
	*/
	if (strncmp(line,"</listOfSpecies>",16) == 0) {
	  in_list_of_species = 0;
	}  else {
	  if (strncmp(line,"<species",8) == 0) {
	    in_species_tag = 1;
	    multiplier = 1.0;
	    line +=8; /* Caution address arithmetic */
	    /*
	      Skip over white space.
	    */
	    nb = count_ws(line);
	    line += nb; /* Caution address arithmetic */
	  }
	}
      }
      if (in_species_tag) {
        ns = count_nws(line);
        while (ns != 0) {
          /*
            Need a routine to scan for keyword "=" "string" triples returning
            full length allowing whitespace before keyword, between any of the
            three tokens and replacing all whitespace in string with 
            underscores, sbml_key_value.
          */
          /*
            check for end tag.
          */
          if ((strncmp(line,"/>",2)==0) || (line[0] == '>')) {
            ns = 0;
            in_species_tag = 0;
            /*
	      Generate line for concs.in file.
	      First convert the amount to a concentration.
	      Also for sbml files that have names where the only
	      real chemical species information is stored 
	      (i.e. the YLI_y7.2.xml file) we need to print
	      a (species_id, species_name) file for analysis and
	      assignment of json id's.
            */
            if (comp[0] != '\0') {
	      /*
		We had a non null compartment.
	      */
	      init_amt *= multiplier * recip_volume;
	      fprintf(concs_fp,"%s:%s\t%le\t%c\n",
		      tspecies,comp,init_amt,vc[variable]);
            } else {
	      /*
		Compartment did not exist.
	      */
	      init_amt *= default_recip_volume;
	      fprintf(concs_fp,"%s\t%le\t%c\n",
		      tspecies,init_amt,vc[variable]);
            }
	    if (species_name[0] != '\0') {
	      fprintf(id_name_fp,"%s\t%s\n",species,species_name);
	      species_name[0] = '\0';
	    }
	    multiplier = 1.0;
	    recip_volume = default_recip_volume;
          } else {
            /*
	      Not the end of species  tag.
            */
            tl = sbml_key_value(line,key,value,max_key_len,max_val_len);
            if (tl <= 0) {
              ns = 0;
            } else {
              line += tl;
              tag = sbml_lookup_species_attribute(key);
              if (tag >= 0) {
		switch (tag) {
		case 0:
		  /*
		    boundary_condition: not used yet.
		  */
		  strcpy (bc,value);
		  break;
		case 1:
		  /*
		    compartment:
		  */
		  strcpy(comp,value);
		  ci = compartment_lookup(comp,&pseudo_state);
		  if (ci >= 0) {
		    compartment = &state->sorted_compartments[ci];
		    recip_volume = compartment->recip_volume;
		  } else {
		    recip_volume = default_recip_volume;
		  }
		  break;
		case 2:
		  /*
		    constant
		  */
		  if (strcmp(value,"true") == 0) {
		    variable = 0;
		  } else {
		    variable = 1;
		  }
		  break;
		case 3:
		  /*
		    hasOnlySubstanceUnits
		    We do not actually use this now,
		    we assume this is false.
		  */
		  if (strcmp(value,"true") == 0) {
		    multiplier = state->recip_avogadro;
		  }
		  break;
		case 4:
		  /*
		    Species id
		  */
		  strcpy(species,value);
		  if (translate_from_modelseed) {
		    tspecies = translate_ms_2_js(ms2js_data,species);
		    if (tspecies == NULL) {
		      if (lfp) {
			fprintf(lfp,"sbml_process_list_of_species: Warning:"
				" species %s did not have a json translation"
				" using as is.\n",
				species);
			fflush(lfp);
		      }
		      tspecies = species;
		    }
		  } else {
		    tspecies = species;
		  }

		  break;
		case 5:
		  /*
		    initialAmount
		  */
		  sscanf(value,"%le",&init_amt);
		  break;
		case 6:
		  /*
		    Species name
		    Here this is not part of the sbml standard
		    but the xml file we got this was the only
		    place the actual chemical name was mentioned.
		    The 
		  */
		  strcpy(species_name,value);
		  break;
		case 7:
		  /*
		    substanceUnits
		  */
		  strcpy(units,value);
		  if (strcmp(value,"mole") == 0)  {
		    multiplier = mole;
		  } else {
		    if (strcmp(value,"millimole") == 0) {
		      multiplier = millimole;
		    } else {
		      if (strcmp(value,"micromole") == 0) {
			multiplier = micromole;
		      } else {
			if (strcmp(value,"nanomole") == 0) {
			  multiplier = nanomole;
			} else {
			  if (strcmp(value,"picomole") == 0) {
			    multiplier = picomole;
			  } else {
			    if (strcmp(value,"femtomole") == 0) {
			      multiplier = femtomole;
			    } else {
			      fprintf(error_fp,"sbml_process_list_of_species: "
				      "Error, unrecognized unit for ammount "
				      "was %s\n",value);
			      fflush(error_fp);
			      multiplier = mole;
			    }
			  }
			}
		      }
		    }
		  }
		  break;
		} /* end switch */  
	      } else {
		fprintf(error_fp,"sbml_process_list_of_species: "
			"Error found unexpected key :%s\n",key);
		fflush(error_fp);
	      }
            } /* end else we found a keyword=value triple. */
          } /* end else we did not find and end species tag */
        } /* end while (ns != 0) */
      } /* end if (in_species_tag) */
    } /* end while (!feof...) */
  } /* end if open of concs.in file succeeded. */
  return (success);
}
