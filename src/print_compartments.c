/* print_compartments.c
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

#include "print_compartments.h"
int  print_compartments(struct state_struct *state) {
  /* 
    print the compartment volumes, ph and ionic_strengths
    Return 1 on success, 0 if the cmpts_echo_file could not be opened
    for writing.

    Called by echo_inputs

    Arguments:
    
    Name          TMF      Description

    state         G*I      state structure :
                           input fields are cmpts_echo_file,
			                    nunique_compartments,
			                    sorted_compartments,
					    compartment_text,
					    lfp.
                           no fields of state are modified.
    
  */
  struct compartment_struct *compartment;
  double volume;
  double ph;
  double ionic_strength;

  char *cmpts_echo_file;
  char *name;
  char *compartment_text;

  int64_t string;
  
  int nunique_compartments;
  int i;

  int success;
  int padi;

  FILE *cmpts_echo_fp;
  FILE *lfp;

  success              = 1;
  lfp                  = state->lfp;
  cmpts_echo_file      = state->cmpts_echo_file;
  cmpts_echo_fp        = fopen(cmpts_echo_file,"w");
  if (cmpts_echo_fp == NULL) {
    success = 0;
    if (lfp) {
      fprintf(lfp,"print_compartments: unable to open %s\n",cmpts_echo_file);
      fflush(lfp);
    }
  }
  if (success) {
    nunique_compartments = state->nunique_compartments;
    compartment          = state->sorted_compartments;
    compartment_text     = state->compartment_text;
    /*
      First print and header and then 
      the default compartment values.
    */
    volume               = compartment->volume;
    ph                   = compartment->ph;
    ionic_strength       = compartment->ionic_strength;
    fprintf(cmpts_echo_fp,"Compartment\tVolume\tpH\tIonic_strength\n"
	    "default\t%le\t%le\t%le\n",volume,ph,ionic_strength);
  
    for (i=1;i<nunique_compartments;i++) {
      string = compartment->string;
      name   = (char *)&compartment_text[string];
      fprintf(cmpts_echo_fp,"%s\t%le\t%le\t%le\n",
	      name,volume,ph,ionic_strength);
      compartment += 1; /* Caution address arithmetic here. */
    }
    fclose(cmpts_echo_fp);
  }
  return(success);
}
