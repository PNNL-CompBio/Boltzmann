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
int compute_molecule_dg0tfs(struct formation_energy_struct *fes){
/*
  Compute the dg0tf for all of the molecules in the reactions file.
  Called by: boltzmann_init.
  Calls:     compute_molecule_dg0tf, fopen, fprintf, fclose, fflush
*/
  struct pseudoisomer_struct *pseudoisomers;
  int num_cpds;

  struct molecule_struct *cur_molecules;
  char *pseudoisomer_strings;
  char *molecules_text;
  char *molecule;

  int success;
  int i;
  int nu_molecules;
  int print_output;

  double ph;
  double temp_kelvin;
  double ionic_strength;
  double m_rt;
  double m_r_rt;
  double min_molecule_dg0tf;
  double *molecule_dg0tfs;


  FILE *lfp;
  success = 1;
  lfp = fes->log_fp;
  pseudoisomers            = fes->pseudoisomers;
  pseudoisomer_strings     = fes->pseudoisomer_strings;
  num_cpds                 = (int)fes->num_pseudoisomers;
  
  nu_molecules             = (int)fes->nunique_molecules;
  cur_molecules            = fes->sorted_molecules;
  molecules_text           = fes->molecules_text;

  molecule_dg0tfs          = fes->molecule_dg0tfs;

  ph                       = fes->ph;
  temp_kelvin              = fes->temp_kelvin;
  ionic_strength           = fes->ionic_strength;
  m_rt                     = fes->m_rt;
  m_r_rt                   = fes->m_r_rt;
  print_output             = (int)fes->print_output;
  if (print_output) {
    if (lfp) {
      fprintf(lfp,"\nOutput from compute_molecule_dg0tfs.c:\n\n");
      fprintf(lfp, "index of molecule in sorted_molecules \t"
	      "molecule name\t deltaG0_tf(kJ/mol)\n");
    }
  }
  for (i=0;i<nu_molecules;i++) {

    molecule    = (char *)&molecules_text[cur_molecules->string];

    compute_molecule_dg0tf(ph,
			   m_rt,
			   m_r_rt,
			   ionic_strength,
			   molecule,
			   pseudoisomers,
			   pseudoisomer_strings,
			   num_cpds,
			   &molecule_dg0tfs[i]);
    if ((i == 0) || (molecule_dg0tfs[i] < min_molecule_dg0tf)) {
      min_molecule_dg0tf = molecule_dg0tfs[i];
    }
    if (print_output) {
      if (lfp) {
	fprintf(lfp,"%i \t %s \t %f \n",i,molecule,
		molecule_dg0tfs[i]);
      }
    }
    cur_molecules += 1; /* Caution address arithmetic. */
  }
  fes->min_molecule_dg0tf = min_molecule_dg0tf;
  if (print_output) {
    if (lfp) {
      fprintf(lfp,"\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n");
    }
  }
  return (success);
}
