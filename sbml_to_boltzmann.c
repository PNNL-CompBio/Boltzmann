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
#include "sbml2bo_structs.h"

#include "sbml_alloc0.h"
#include "sbml_set_file_names.h"
#include "sbml_alloc1.h"
#include "sbml_count_cmpts.h"
#include "parse_sbml.h"

#include "sbml_to_boltzmann.h"

int sbml_to_boltzmann(struct state_struct *state) {
  struct sbml2bo_struct *sbml_state;
  struct compartment_struct *unsorted_compartments;
  struct compartment_struct *sorted_compartments;
  int success;
  int max_compartment_len;
  FILE *lfp;
  FILE *extra_fp;
  success = 1;
  /*
    Need to allocate space for the state structure and its
    filename fields.
  */
  success = sbml_alloc0(&sbml_state);
  /*
    Need to read input file name and generate output file names.
  */
  if (success) {
    strcpy(sbml_state->sbml_file,state->sbml_file);
    success = sbml_set_file_names(sbml_state);
  }
  /*
    Now we need to actually count the compartment tabs in the sbml file to
    get a handle on the number of compartments and allocate
    an array to store comparment names and volumes in order to sort them
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
    strcpy(state->compartment_file,sbml_state->cmpts_dat_file);          
  }
  return(success);
}  
