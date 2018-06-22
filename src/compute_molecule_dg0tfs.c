/* compute_molecule_dg0tfs.c
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
 * compute_molecule_dg0tf.c
 *
 *  Created on: Apr 18, 2013
 *      Author:  Dennis G. Thomas
        Modified by Doug Baxter June 2, 2013
*/
#include "boltzmann_structs.h"

#include "compute_molecule_dg0tf.h"

#include "compute_molecule_dg0tfs.h"
int compute_molecule_dg0tfs(struct state_struct *state,
			    struct pseudoisomer_struct *pseudoisomers,
			    char *pseudoisomer_strings,
			    int  num_cpds) {
/*
  Compute the dg0tf for all of the molecules in the reactions file.
  num_cpds is the number of entries in the pseudoismers array.
  Called by: compute_standard_energies
  Calls:     compute_molecule_dg0tf, fopen, fprintf, fclose, fflush
*/

  struct molecule_struct *cur_molecules;
  struct compartment_struct *compartments;
  struct compartment_struct *compartment;
  char *molecules_text;
  char *molecule_name;
  char *compartment_text;
  char *compartment_name;

  double ph;
  double temp_kelvin;
  double ionic_strength;
  double m_rt;
  double m_r_rt;
  double min_molecule_dg0tf;
  double *molecule_dg0tfs;

  int *dg0tfs_set;
  int success;
  int print_output;

  int i;
  int nu_molecules;

  int use_dgzero;
  int c_index;

  int found;
  int use_pseudoisomers;

  FILE *lfp;
  FILE *efp;

  success            = 1;
  lfp                = state->lfp;
  use_dgzero         = state->use_dgzero;
  use_pseudoisomers  = state->use_pseudoisomers;
  dg0tfs_set         = state->dg0tfs_set;
  nu_molecules       = (int)state->nunique_molecules;
  cur_molecules      = state->sorted_molecules;
  compartments       = state->sorted_compartments;
  molecules_text     = state->molecules_text;
  compartment_text   = state->compartment_text;

  molecule_dg0tfs    = state->molecule_dg0tfs;
  /*
  ph                 = state->ph;
  ionic_strength     = state->ionic_strength;
  */
  temp_kelvin        = state->temp_kelvin;
  m_rt               = state->m_rt;
  m_r_rt             = state->m_r_rt;
  print_output       = (int)state->print_output;
  min_molecule_dg0tf = 0.0;
  if (print_output) {
    if (lfp) {
      fprintf(lfp,"\nOutput from compute_molecule_dg0tfs.c:\n\n");
      fprintf(lfp, "index of molecule in sorted_molecules \t"
	      "molecule name:compartment\t deltaG0_tf(kJ/mol)\n");
    }
  }
  for (i=0;((i<nu_molecules) && success);i++) {

    molecule_name    = (char *)&molecules_text[cur_molecules->string];
    c_index          = cur_molecules->c_index;
    compartment      = (struct compartment_struct *)&compartments[c_index];
    compartment_name = (char *)&compartment_text[compartment->string];

    found = 0;
    if (use_pseudoisomers == 2) {
      if (cur_molecules->potential_set) {
	molecule_dg0tfs[i] = cur_molecules->potential;
	found = 1;
      }
    }
    if (found == 0) {
      ionic_strength   = compartment->ionic_strength;
      ph               = compartment->ph;
      found = compute_molecule_dg0tf(ph,
				     m_rt,
				     m_r_rt,
				     ionic_strength,
				     molecule_name,
				     pseudoisomers,
				     pseudoisomer_strings,
				     num_cpds,
				     &molecule_dg0tfs[i]);
    }
    if (found) {
      dg0tfs_set[i] = 1;
      if ((i == 0) || (molecule_dg0tfs[i] < min_molecule_dg0tf)) {
	min_molecule_dg0tf = molecule_dg0tfs[i];
      }
      if (print_output) {
	if (lfp) {
	  if (c_index > 0) {
	    fprintf(lfp,"%i \t %s:%s \t %f \n",i,molecule_name,
		    compartment_name,molecule_dg0tfs[i]);
	  } else {
	    fprintf(lfp,"%i \t %s \t %f \n",i,molecule_name,
		    molecule_dg0tfs[i]);
	  }
	}
      }
    } else {
      dg0tfs_set[i] = 0;
      /*
	Check to see if use_dgzero is set, if so then reactions involving
	this molecule will use the read in DGZERO value from the reaction.dat 
	file. Otherwise it is an error, set success to 0 and print messsage
	to log file and quit.
      */
      if (use_dgzero) {
	if (lfp) {
	  fprintf(lfp,"compute_molcule_dg0tfs: %s not found in pseudoisomer file, using specified DGZERO from reactions.dat file for reactions it is involved in\n",molecule_name);
	  fflush(lfp);
	}
      } else {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"compute_molecule_dg0tfs: %s not found in pseudoisomer file and use_dgzero not specified. Quitting.\n",molecule_name);
	  fflush(lfp);
	}
	break;
      }
    }
    cur_molecules += 1; /* Caution address arithmetic. */
  }
  state->min_molecule_dg0tf = min_molecule_dg0tf;
  if (print_output) {
    if (lfp) {
      fprintf(lfp,"\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n");
    }
  }
  return (success);
}
