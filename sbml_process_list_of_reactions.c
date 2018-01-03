/* sbml_process_list_of_reactions.c
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

#include "sbml2bo_structs.h"

#include "count_ws.h"
#include "count_nws.h"
#include "sbml_look_for_in_reaction_tag.h"
#include "sbml_key_value.h"
#include "sbml_lookup_reaction_attribute.h"
#include "sbml_lookup_speciesref_attribute.h"

#include "sbml_process_list_of_reactions.h"

int sbml_process_list_of_reactions(FILE *sbml_fp,
				   char *sbml_buffer,
				   int sbml_buffer_len, 
				   struct sbml2bo_struct *state) {
  /*
    Process the <listOfReactions> tag contents, and generate
    the rxns.dat file.
    Called by: parse_sbml
    Calls:     fprintf, feof, fgets, fflush, strcpy, strcmp, strncmp, sscanf
               count_ws, count_nws,
	       sbml_look_for_in_reaction_tag,
	       sbml_key_value,
               sbml_lookup_reaction_attribute
               sbml_lookup_speciesref_attribute
	       
  */
  /*
    Thoughts:
    We have reaction groups, headed by a rxn_header tag, 
    "<reaction" {in_reaction=1,in_reaction_tag=1}key_value_triples (id, reversible, fast, compartment) ">" {in_reaction_tag=0, generate REACTION line and compartment line, comp[0] = '\0', rxn_count += 1}
    "<listOfReactants" {in_list_of_reactants=1,in_lor_tag=1,species_count=0} key_value_triple (compartment)">" {in_lor_tag = 0,if compartment print LEFT_COMPARTMENT comp, comp[0]='\0',print LEFT, sc = 0}
    "<speciesReference" {in_species_ref_tag = 1} key_value_triple  
                           (species, stoichiometry, constant, compartment) "/>"
			{if (sc == 0) print '+'; print stoic species:compartment		            ;sc += 1. in_species_ref_tag = 0, comp[0]='\0'}

  "</listOfReactants>" {in_list_of_reactants = 0, sc = 0
    "<listOfProducts" {in_list_of_products=1, in_lop_tag=1 sc = 0} key_value_triple (compartment) ">" {in_lop_tag = 0, if compartment print RIGHT_COMPARTMENT comp,comp[0] = '\0', print RIGHT, sc = 0}
    "</listOfProducts>" {in_list_of_prodeuct = 0, sc= 0}
    "<kineticLaw" key_value_triple ">" ignore lines from here		
    "</kineticLaw>" to here
    "</reaction>" {in_reaction=0}
    {print //}
  */
  char comp_c[1024];
  char species_c[1024];
  char value_c[1024];
  char rxn_id_c[1024];
  char key_c[1024];
  char units_c[1024];
  char *key;
  char *value;
  char *comp;
  char *species;
  char *rxn_id;
  char *units;
  char bcf[8];
  char vc[8];
  char *concs_file;
  char *bc;
  char *line;
  
  double init_amt;
  double delta_g0;

  int variable;
  int in_reaction_tag;

  int in_list_of_reactions;
  int in_list_of_reactants;

  int in_list_of_products;
  int in_kinetic_law;
  
  int success;
  int rxn_count;

  int enclosing_tag;
  int not_in_tag;

  int in_list_of_reactants_tag;
  int end_list_of_reactants_tag;

  int in_list_of_products_tag;
  int end_list_of_products_tag;

  int in_species_reference_tag;
  int coefficient;

  int in_kinetic_law_tag;
  int end_reaction_tag;

  int in_reaction;
  int species_count;

  int nb;
  int ns;

  int tl;
  int tag;

  int max_key_len;
  int max_val_len;

  FILE *lfp;
  FILE *error_fp;
  FILE *rxns_fp;

  success     = 1;
  key         = (char *)&(key_c[0]);
  value       = (char *)&(value_c[0]);
  comp        = (char *)&(comp_c[0]);
  species     = (char *)&(species_c[0]);
  units       = (char *)&(units_c[0]);
  rxn_id      = (char *)&(rxn_id_c[0]);
  bc          = (char *)&bcf[0];
  vc[0]       = 'F';
  vc[1]       = 'V';
  max_key_len = 1023;
  max_val_len = 1023;

  lfp        = state->log_fp;
  rxns_fp    = state->rxns_fp;
  if (lfp == NULL) {
    error_fp = stderr;
  } else {
    error_fp = lfp;
  }
  /*
    Set the default delta_g0 and delta_g0_units.
  */
  delta_g0 = 0.0;
  strcpy(units,"KJ/MOL");

  comp[0] = '\0';
  species[0] = '\0';
  init_amt=0;
  units[0] = '\0';
  variable     = 1;
  in_reaction = 0;
  in_reaction_tag = 0;
  rxn_count          = 0;
  not_in_tag                = 1;
  in_reaction_tag           = 2;
  in_list_of_reactants_tag  = 4;
  end_list_of_reactants_tag = 8;
  in_list_of_products_tag   = 16;
  end_list_of_products_tag  = 32;
  in_species_reference_tag  = 64;
  in_kinetic_law_tag        = 128;
  end_reaction_tag          = 256;

  enclosing_tag = not_in_tag;
  species_count = 0;
  while (!feof(sbml_fp) && (in_list_of_reactions == 1)) {
    /* 
       Now we are in the listOfReactions section.    
    */
    fgets(sbml_buffer,sbml_buffer_len,sbml_fp);
    line = sbml_buffer;
    /*
      Strip off any preceding white space.
    */
    nb = count_ws(line);
    line += nb; /* Caution address arithmetic */
    /*
      Check if we are not in a reaction.
    */
    if (in_reaction == 0) {
      /*
	First check for an </listOfReactions> tag.
      */
      if (strncmp(line,"</listOfReactions>",18) == 0) {
	in_list_of_reactions = 0;
      }  else {
	/*
	  look for a <reaction start to the line.
	*/
	if (strncmp(line,"<reaction",9) == 0) {
	  in_reaction = 1;
	  enclosing_tag = in_reaction_tag;
	  line +=9; /* Caution address arithmetic */
	}
      }
    }
    if (in_reaction) {
      if (enclosing_tag == not_in_tag) {
	enclosing_tag = sbml_look_for_in_reaction_tag(line);
      }
      switch (enclosing_tag) {
      case 1 :
	/*
	  We are not currently in a tag, skip line
	  we might want to log skipped lines.
	*/
	break;
      case 2:
	/*
	  in_reaction_tag
	  "<reaction" tab has already been skipped over. 
	  Look for end of reaction tag,">"
	*/
	nb = count_ws(line);
	line += nb; /* Caution address arithmetic */
	ns = count_nws(line);
	while (ns != 0) {
	  if (line[0] == '>') {
	    /*
	      We have hit the end of the reaction tag,
	      generate  the reaction label, and
	      compartment if there was one, reset
	      the compartment, enclosing_tag and set the
	      species count to 0. and break out of this loop.
	    */
	    fprintf(rxns_fp,"REACTION\t%s\n",rxn_id);
	    if (comp[0] != '\0') {
	      fprintf(rxns_fp,"COMPARTEMNT\t%s\n",comp);
	    }
	    comp[0] = '\0';
	    species_count = 0;
	    enclosing_tag = not_in_tag;
	    ns = 0;
	  } else {
	    /*
	      look for key=value triples, expecting keys of
	      "id", "reversible" "fast" and possibly "compartment"
	    */
	    tl = sbml_key_value(line,key,value,max_key_len,max_val_len);
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
		} /* end switch */  
	      } else {
		fprintf(error_fp,"sbml_process_list_of_reactions: "
			"Error found unexpected reaction key :%s\n",key);
		fflush(error_fp);
	      }
	    } /* end else we found a keyword=value triple. */
	  } /* end else we did not find and en reaction_tag tag */
	} /* end while (ns != 0), in reaction tag. */
	break;
      case 4:
	/*
	  In list_of_reactants tag.
	  skip over "<listOfReactants"
	*/
	line += 16;
	nb = count_ws(line);
	line += nb; /* Caution address arithmetic */
	ns = strlen(line);
	while (ns != 0)  {
	  if (line[0] == '>') {
	    /*
	      End of list_of_reactants tab,
	      print a LEFT_COMPARTMENT line if compartment is not null,
	      Then print the start of the LEFT line.
	    */
	    if (comp[0] != '\0') {
	      fprintf(rxns_fp,"LEFT_COMPARTMENT\t%s\n",comp);
	      comp[0] = '\0';
	    }
	    fprintf(rxns_fp,"LEFT\t");
	    enclosing_tag = not_in_tag;
	    species_count = 0;
	    ns = 0;
	  } else {
	    tl = sbml_key_value(line,key,value,max_key_len,max_val_len);
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
		fprintf(error_fp,"sbml_process_list_of_reactions: Error "
			"Unrecognized key=value pair in listOfReactants "
			"tag:%s=%s\n",key,value);
	      }
	    } /* end else found a key=value field in the listOfReactans tag.*/
	  } /* end else didn't find end of listOfReactants tag */	  
	} /* end while (ns != 0) */
	break;
      case 8:
	/*
	  We found the </listOfReactants> tag, write a new line and
	  set the species count to zero and the compartment to the 
	  null string.
	*/
	fprintf(rxns_fp,"\n");
	species_count = 0;
	comp[0] = '\0';
	enclosing_tag = not_in_tag;
	break;
      case 16:
	/* 
	   In list_of_products tag.
	   We found the <listOfProducts> tag look for a compartment.
	   skip over "<listOfProducts"
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
	    tl = sbml_key_value(line,key,value,max_key_len,max_val_len);
	    if (tl <= 0) {
	      ns = 0;
	    } else {
	      line += tl;
	      nb = count_ws(line);
	      line += nb; /* Caution address arithmetic */
	      ns = strlen(line);
	      if (strcmp(key,"compartment") == 0) {
		strcpy(comp,value);
	      } else {
		fprintf(error_fp,"sbml_process_list_of_reactions: Error "
			"Unrecognized key=value pair in listOfProducts "
			"tag:%s=%s\n",key,value);
	      }
	    } /* end else found a key=value field in the listOfProducts tag.*/
	  } /* end else didn't find end of listOfProducts tag */	  
	} /* end while (ns != 0) */
	break;
      case 32:
	/*
	  We found the </listOfProducts> tag, write a new line and
	  set the species count to zero and the compartment to the 
	  null string.
	*/
	fprintf(rxns_fp,"\n");
	species_count = 0;
	comp[0] = '\0';
	enclosing_tag = not_in_tag;
	break;
      case 64:
	/*
	  We found a "<speciesReference" tag, process the fields their
	  storing the species name, its stoichiometric coefficient,
	  its compartment if specified and a constant field ignored
	  here.
	*/
	line += 17;
	nb = count_ws(line);
	line += nb; /* Caution address arithmetic */
	ns = count_nws(line);
	while (ns != 0) {
	  if (strncmp(line,"/>",2) == 0) {
	    /*
	      We have reached the end of the speciesReference tag,
	      print out the stoichiometry, species name and compartment
	      if it exists, and preceding + sign if species_count is > 0.
	    */
	    if (species_count > 0) {
	      fprintf(rxns_fp,"\t+\t");
	    }
	    if (coefficient > 1) {
	      fprintf(rxns_fp,"%d\t",coefficient);
	    }
	    fprintf(rxns_fp,"%s",species);
	    if (comp[0] != '\0') {
	      fprintf(rxns_fp,":%s",comp);
	      comp[0] = '\0';
	    }
	    species_count += 1;
	    enclosing_tag = not_in_tag;
	    ns = 0;
	  } else {
	    /*
	      look for key=value triples, expecting keys of
	      "species", "stoichiometry" "constant" and possibly "compartment"
	    */
	    tl = sbml_key_value(line,key,value,max_key_len,max_val_len);
	    if (tl <= 0) {
	      ns = 0;
	    } else {
	      line += tl;
	      nb = count_ws(line);
	      line += nb; /* Caution address arithmetic */
	      ns = strlen(line);
	      line += tl;
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
		  break;
		case 3:
		  /*
		    stoichiomentry.
		  */
		  sscanf(value,"%d",&coefficient);
		  break;
		} /* end switch */  
	      } else {
		fprintf(error_fp,"sbml_process_list_of_reactions: "
			"Error found unexpected speciesref key :%s\n",key);
		fflush(error_fp);
	      } 
	    } /* end else we had key=value triple */
	  } /* end else not at end of speciesReference tag. */
	} /* end while (ns ...) */
	break;
      case 128:
	/*
	  in kineticLaw tag, ignore lines till we reach the end.
	*/
	nb = count_ws(line);
	line += nb; /* Caution address arithmetic */
	ns = count_nws(line);
	if (strncmp(line,"</kineticLaw>",13) == 0) {
	  enclosing_tag = 0;
	}
	break;
      case 256:
	enclosing_tag = 0;
	in_reaction = 0;
	fprintf(rxns_fp,"DGZERO\t%le\n",delta_g0);
	fprintf(rxns_fp,"DGZERO_UNITS\t%s\n",units);
	/*
	  Reset the default delta_g0 and delta_g0_units.
	*/
	fprintf(rxns_fp,"//\n");
	delta_g0 = 0.0;
	strcpy(units,"KJ/MOL");
	rxn_count += 1;
      } /* end switch(enclosing_tag) */
    } /* end if_in_reaction */
  } /* end while (!feof...) */
  state->num_reactions = rxn_count;
  return (success);
}
