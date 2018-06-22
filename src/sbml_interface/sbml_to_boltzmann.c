/* sbml_to_boltzmann.c
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
/*
  Need to extract information from two maybe three sections.
  Want to produce a concs.in rxns.dat and maybe a compartments.dat 
  file for input to boltzmann.
  Will need the name of an sbml input file,
  the names of the three output files, Possibly some more parameters.
  Might also want to have a output .lng file with the number of compartments,
  number of species, and number of reactions in it.

  main program might well be called sbml2bo short for sbml to boltzmann

  I think we need a state struct, let us say sbml2bo_struct, to
  hold input parameters, file pointers, counters and general state.

  List of compartments, notion of a volume per compartment
  we don't have that yet but maybe ought to. Ask Bill about this.
  tags of interest: might need listOfUnitDefinitions but not at first.
    id
    size (maybe)
    spatialDimensions (maybe)
    units(maybe)
    constant (maybe)
  If we don't want volumes we'll get compartments from reaction list.  

  The listOfSpecies section to build tine initial concentrations file.
      tags of interest
          comp,
          id,
	  intialAmount
          substanceUnits (maybe)
	  constant
   
  May want a version that has substance units.
  The listOfReactions section.
    tags of interest
     reaction id
     listOfReactants
             speciesReference
	          Species
		  stoichiometry
		  constant
     listOfProducts 
             speciesReference
	          Species
		  stoichiometry
		  constant

      these two unlikeliy:
     KineticLaw if it has a k_eq
     listOfLocalParameters if they have a dg0
*/
#include "boltzmann_structs.h"

#include "size_ms2js_file.h"
#include "size_kg2js_file.h"
#include "sbml_alloc0.h"
#include "sbml_set_file_names.h"
#include "sbml_alloc2.h"
#include "read_ms2js.h"
#include "read_kg2js.h"
#include "sort_json_ids.h"
#include "sbml_count_cmpts.h"
#include "sbml_count_species.h"
#include "sbml_alloc1.h"
#include "parse_sbml.h"

#include "sbml_to_boltzmann.h"

int sbml_to_boltzmann(struct state_struct *state) {
  /*
    Convert an sbml input file to a rxns.dat and concs.in files.
    Called by: io_size_init
    Calls:     size_ms2js_file.h,
    	       size_kg2js_file.h,
    	       sbml_alloc0.h,
    	       sbml_set_file_names.h,
    	       sbml_alloc2.h,
    	       read_ms2js.h,
    	       read_kg2js.h,
    	       sort_json_ids.h,
    	       sbml_count_cmpts.h,
    	       sbml_count_species.h,
    	       sbml_alloc1.h,
    	       parse_sbml.h
  */
  struct sbml2bo_struct *sbml_state;
  int64_t num_modelseed_ids;
  int64_t num_kegg_ids;
  int64_t length_ms2js_strings;
  int64_t length_kg2js_strings;
  int success;
  int padi;
  FILE *lfp;
  FILE *extra_fp;
  success = size_ms2js_file(state,&num_modelseed_ids,
			    &length_ms2js_strings);
  
  if (success) {
    success = size_kg2js_file(state,&num_kegg_ids,
			      &length_kg2js_strings);
  }
  /*
    Need to allocate space for the state structure and its
    filename fields and the ms2js structure and its fields,
    ms2js_strings, ms_ids, js_ids.
  */
  if (success) {
    success = sbml_alloc0(&sbml_state);
  }
  if (success) {
    strcpy(sbml_state->sbml_file,state->sbml_file);
    /*
      Need to generate output file names.
    */
    success = sbml_set_file_names(sbml_state);
  }
  if (success) {
    strcpy(sbml_state->ms2js_file,state->ms2js_file);
    strcpy(sbml_state->kg2js_file,state->kg2js_file);
    sbml_state->num_modelseed_ids    = num_modelseed_ids;
    sbml_state->num_kegg_ids         = num_kegg_ids;
    sbml_state->length_ms2js_strings = length_ms2js_strings;
    sbml_state->length_kg2js_strings = length_kg2js_strings;
    sbml_state->default_comp_size    = 1.0e-15;
    success = sbml_alloc2(sbml_state,num_modelseed_ids,length_ms2js_strings,
			  num_kegg_ids,length_kg2js_strings);
  }
  if (success) {
    success = read_ms2js(sbml_state);
  }
  if (success) {
    success = read_kg2js(sbml_state);
  }
  if (success) {
    /*
      We also need to have a sorted list of json_id's to look up.
      We will use those corresponding to the kegg_ids as the list.
    */
    success = sort_json_ids(sbml_state);
  }
  /*
    Now we need to actually count the compartment tabs in the sbml file to
    get a handle on the number of compartments and allocate
    an array to store compartment names and volumes in order to sort them
    and be able to look them up when processing the list_of_species tab
    in computing concentrations from the amounts and corresponding 
    compartment volumes. We will assume that the first listed compartment
    in the listOfCompartments is the default compartment when not specified.
  */
  if (success) {
    /*
      Count the compartments, setting the num_cmpts field in sbml_state.
    */
    success = sbml_count_cmpts(sbml_state);
  }
  if (success) {
    success = sbml_count_species(sbml_state);
  }
  if (success) {
    success = sbml_alloc1(sbml_state);
  }
  if (success) {
    /*
      Set log file from state.
    */
    lfp = state->lfp;
    sbml_state->log_fp = lfp;
    success = parse_sbml(sbml_state);
  }
  if (success) {
    strcpy(state->init_conc_file,sbml_state->concs_in_file);
    strcpy(state->reaction_file,sbml_state->rxns_dat_file);
    strcpy(state->compartments_file,sbml_state->cmpts_dat_file);          
  }
  return(success);
}  
