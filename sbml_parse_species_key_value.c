/* sbml_parse_species_key_value.c
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
#include "sbml_lookup_species_attribute.h"
#include "compartment_lookup.h"
#include "sbml_process_substanceunits.h"

#include "sbml_parse_species_key_value.h"
int sbml_parse_species_key_value(struct sbml2bo_struct *state,
				 char *line,
				 double default_recip_volume,
				 char *bc,
				 char *comp,
				 char *species,
				 char *species_name,
				 double *recip_volume,
				 double *multiplier,
				 double *init_amt,
				 double *init_conc,
				 int *variable) {
  /*
    Interpret the attribute value of a key_value attribute pair
    in a <species> tag. Returns the total length of the key and value
    pair in the line buffer.
    Called by: sbml_process_list_of_species.
    Calls:     sbml_read_key_value;            
               sbml_lookup_species_attribute,
               compartment_lookup,
	       sbml_process_substanceunits,
	       strcpy, strcmp, fprintf, fflush

  */
  struct state_struct pseudo_state;
  struct compartment_struct *compartment;
  char key_c[1024];
  char value_c[1024];
  char units_c[1024];
  char *key;
  char *value;
  char *units;

  int max_key_len;
  int max_val_len;

  int tag;
  int ci;

  int tl;
  int padi;

  FILE *lfp;
  FILE *error_fp;
  key         = (char *)&(key_c[0]);
  value       = (char *)&(value_c[0]);
  units       = (char *)&(units_c[0]);
  lfp = state->log_fp;
  if (lfp == NULL) {
    error_fp = stderr;
  } else {
    error_fp = lfp;
  }
  max_key_len = 1023;
  max_val_len = 1023;
  pseudo_state.compartment_text     = state->compartment_text;
  pseudo_state.sorted_compartments  = state->sorted_compartments;
  pseudo_state.nunique_compartments = state->num_cmpts;

  tl = sbml_read_key_value(line,key,value,max_key_len,max_val_len);
  if (tl > 0) {
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
	  *recip_volume = compartment->recip_volume;
	} else {
	  *recip_volume = default_recip_volume;
	}
	break;
      case 2:
	/*
	  constant
	*/
	if (strcmp(value,"true") == 0) {
	  *variable = 0;
	} else {
	  *variable = 1;
	}
	break;
      case 3:
	/*
	  hasOnlySubstanceUnits
	  We do not actually use this now,
	  we assume this is false.
	*/
	if (strcmp(value,"true") == 0) {
	  *multiplier = state->recip_avogadro;
	}
	break;
      case 4:
	/*
	  Species id
	*/
	strcpy(species,value);
	break;
      case 5:
	/*
	  initialAmount
	*/
	sscanf(value,"%le",init_amt);
	break;
      case 6:
	/*
	  initialConcentration
	*/
	sscanf(value,"%le",init_conc);
	break;
      case 7:
	/*
	  Species name
	  Here this is not part of the sbml standard
	  but the xml file we got this was the only
	  place the actual chemical name was mentioned.
	  The 
	*/
	strcpy(species_name,value);
	break;
      case 8:
	/*
	  substanceUnits
	*/
	strcpy(units,value);
	sbml_process_substanceunits(units,multiplier,error_fp);
	break;
      } /* end switch */  
    } else {
      /*
	Too much information if we print all these out.
      fprintf(error_fp,"sbml_parse_species_key_value: "
	      "Error found unexpected key:%s in species tag\n",key);
      fflush(error_fp);
      */
    }
  }
  return(tl);
}
