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
  char *molecules_text;
  char *molecule;

  double ph;
  double temp_kelvin;
  double ionic_strength;
  double m_rt;
  double m_r_rt;
  double min_molecule_dg0tf;
  double *molecule_dg0tfs;

  int success;
  int print_output;

  int i;
  int nu_molecules;

  FILE *lfp;
  FILE *efp;

  success = 1;
  lfp = state->lfp;
  
  nu_molecules             = (int)state->nunique_molecules;
  cur_molecules            = state->sorted_molecules;
  molecules_text           = state->molecules_text;

  molecule_dg0tfs          = state->molecule_dg0tfs;

  ph                       = state->ph;
  temp_kelvin              = state->temp_kelvin;
  ionic_strength           = state->ionic_strength;
  m_rt                     = state->m_rt;
  m_r_rt                   = state->m_r_rt;
  print_output             = (int)state->print_output;
  min_molecule_dg0tf       = 0.0;
  if (print_output) {
    if (lfp) {
      fprintf(lfp,"\nOutput from compute_molecule_dg0tfs.c:\n\n");
      fprintf(lfp, "index of molecule in sorted_molecules \t"
	      "molecule name\t deltaG0_tf(kJ/mol)\n");
    }
  }
  for (i=0;((i<nu_molecules) && success);i++) {

    molecule    = (char *)&molecules_text[cur_molecules->string];

    success = compute_molecule_dg0tf(ph,
				     m_rt,
				     m_r_rt,
				     ionic_strength,
				     molecule,
				     pseudoisomers,
				     pseudoisomer_strings,
				     num_cpds,
				     &molecule_dg0tfs[i]);
    
    if (success) {
      if ((i == 0) || (molecule_dg0tfs[i] < min_molecule_dg0tf)) {
	min_molecule_dg0tf = molecule_dg0tfs[i];
      }
      if (print_output) {
	if (lfp) {
	  fprintf(lfp,"%i \t %s \t %f \n",i,molecule,
		  molecule_dg0tfs[i]);
	}
      }
    } else {
      break;
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
