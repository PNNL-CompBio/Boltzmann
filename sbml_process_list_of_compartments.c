/* sbml_process_list_of_compartments.c
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
#include "sbml_lookup_compartment_attribute.h"
#include "sbml_volume_units_conversion.h"
#include "sort_compartments.h"

#include "sbml_process_list_of_compartments.h"

int sbml_process_list_of_compartments(FILE *sbml_fp,
				      char *sbml_buffer,
				      int sbml_buffer_len, 
				      struct sbml2bo_struct *state) {
  /*
    Process the <listOfCompartments> tag contents, and generate
    the cmpts.dat file.
    Called by: parse_sbml
    Calls:     fprintf, feof, fgets, fflush, strcpy, strcmp, strncmp,
               count_ws, count_nws,
	       sbml_read_key_value,
               sbml_lookup_compartment_attribute
    Note, it works out that we need to build a list of compartments
    and their voluems/reciprocal volumes for use in the 
    sbml_process_list_of_species section as boltzmann needs intial amounts
    in concentrations and sbml only specifies initial amounts in quantities.
    
  */
  struct compartment_struct *unsorted_compartments;
  struct compartment_struct *sorted_compartments;
  struct compartment_struct *compartment;
  char *compartment_text;
  char comp_c[1024];
  char units_c[1024];
  char key_c[1024];
  char value_c[1024];
  char *key;
  char *value;
  char *comp;
  char *units;
  char bcf[8];
  char vc[8];
  char *concs_file;
  char *bc;
  char *line;
  int64_t cmpt_pos;
  int64_t alignment;
  int64_t align_mask;
  
  double size;
  double multiplier;

  int variable;
  int in_cmpt_tag;

  int in_list_of_cmpts;
  int spatial_dim;
  
  int success;
  int nb;

  int ns;
  int tl;

  int max_key_len;
  int max_val_len;

  int tag;
  int n_cmpts;

  int cmpt_len;
  int padi;
  
  FILE *lfp;
  FILE *error_fp;
  FILE *cmpts_fp;

  success     = 1;
  key         = (char *)&(key_c[0]);
  value       = (char *)&(value_c[0]);
  comp        = (char *)&(comp_c[0]);
  units       = (char *)&(units_c[0]);
  bc          = (char *)&bcf[0];
  max_key_len = 1023;
  max_val_len = 1023;
  vc[0]       = 'F';
  vc[1]       = 'V';
  alignment   = (int64_t)state->alignment;
  align_mask  = (int64_t)state->align_mask;
  cmpt_pos    = (int64_t)0;
  cmpts_fp    = state->cmpts_fp;
  unsorted_compartments = state->unsorted_compartments;
  sorted_compartments = state->sorted_compartments;
  compartment = unsorted_compartments;
  n_cmpts     = 0;
  size        = state->default_comp_size;
  multiplier  = 1.0;
  lfp         = state->log_fp;
  if (lfp == NULL) {
    error_fp = stderr;
  } else {
    error_fp = lfp;
  }
  comp[0] = '\0';
  units[0] = '\0';
  variable     = 1;
  in_cmpt_tag = 0;
  in_list_of_cmpts = 1;
  while (!feof(sbml_fp) && (in_list_of_cmpts == 1)) {
    /* 
  	 Now we are in the listOfCompartments section.    
    */
    fgets(sbml_buffer,sbml_buffer_len,sbml_fp);
    line = sbml_buffer;
    if (in_cmpt_tag == 0) {
  	nb = count_ws(line);
  	line += nb; /* Caution address arithmetic */
  	/*
  	  First check for an </listOfCompartments> tag.
  	*/
  	if (strncmp(line,"</listOfCompartments>",18) == 0) {
  	  in_list_of_cmpts = 0;
  	}  else {
  	  if (strncmp(line,"<compartment",12) == 0) {
  	    in_cmpt_tag = 1;
  	    line +=12; /* Caution address arithmetic */
  	    /*
  	      Skip over white space.
  	    */
  	    nb = count_ws(line);
  	    line += nb; /* Caution address arithmetic */
  	  }
  	}
    }
    if (in_cmpt_tag) {
      ns = count_nws(line);
      while (ns != 0) {
        /*
          Need a routine to scan for keyword "=" "string" triples returning
          full length allowing whitespace before keyword, between any of the
          three tokens and replacing all whitespace in string with 
          underscores, sbml_read_key_value.
        */
        /*
          check for end tag.
        */
        if ((strncmp(line,"/>",2)==0) || (line[0] == '>')) {
          ns = 0;
          in_cmpt_tag = 0;
          /*
  	      Generate line for cmpts.dat file.
          */
	  size = size * multiplier;
	  if (n_cmpts == 0) {
	    state->default_comp_size = size;
	  } 
	  fprintf(cmpts_fp,"%s\t%le\tliters\t%d\t%c\n",
		  comp,size,spatial_dim,vc[variable]);
	  if (n_cmpts < state->num_cmpts) {
	    cmpt_len = strlen(comp) + 1;
	    cmpt_len += ((alignment - (cmpt_len & align_mask)) & align_mask);
	    strcpy((char*)&compartment_text[cmpt_pos],comp);
	    compartment->string = cmpt_pos;
	    compartment->c_index = n_cmpts;
	    compartment->volume = size;
	    cmpt_pos += cmpt_len;
	    n_cmpts += 1;
	    if (size <= 0) {
	      success = 0;
	      ns = 0;
	      fprintf(error_fp,"sbml_process_list_of_compartment: Error, "
		      "0 or negative size for %s\n",comp);
	      fflush(error_fp);
	    } else {
	      compartment->recip_volume = 1.0/size;
	    }
	    compartment += 1; /* Caution address arithmetic here */
	    /*
	      Reset size to default size.
	    */
	  } else {
	    success = 0;
	    ns = 0;
	    fprintf(error_fp,"sbml_process_list_of_compartment: Error, "
		    "Too many compartments. "
		    "Look for bug in sbml_count_cmpts\n");
	    fflush(error_fp);
	  }
        } else {
          /*
  	      Not the end of compartment  tag.
          */
          tl = sbml_read_key_value(line,key,value,max_key_len,max_val_len);
          if (tl <= 0) {
            ns = 0;
          } else {
            line += tl;
            tag = sbml_lookup_compartment_attribute(key);
            if (tag >= 0) {
	      switch (tag) {
	      case 0:
		/*
		  constant
		*/
  	        if (strcmp(value,"true") == 0) {
		  variable = 0;
		} else {
		  variable = 1;
		}
		break;
	      case 1:
		/*
		  id
		*/
		strcpy(comp,value);
		break;
	      case 2:
		/*
		  size
		*/
		sscanf(value,"%le",&size);
		break;
	      case 3:
		/*
		  spatialDimensions
		*/
		sscanf(value,"%d",&spatial_dim);
		break;
	      case 4:
		/*
		  units
		*/
		strcpy(units,value);
		multiplier = sbml_volume_units_conversion(units,error_fp);
		break;
	      } /* end switch */  
	    } else {
	      /*
		This generates too much output, we don't want to know.
	      fprintf(error_fp,"sbml_process_list_of_compartments: "
		      "Error found unexpected key :%s\n",key);
	      fflush(error_fp);
	      */
	    }
          } /* end else we found a keyword=value triple. */
        } /* end else we did not find and end species tag */
      } /* end while (ns != 0) */
    } /* end if (in_compartment_tag) */
  } /* end while (!feof...) */
  /*
    Now sort compartments for use in sbml_process_list_of_species.
  */
  if (success) {
    success = sort_compartments(unsorted_compartments,sorted_compartments,
				compartment_text,n_cmpts);
  }
  state->num_cmpts = n_cmpts;
  return (success);
}
