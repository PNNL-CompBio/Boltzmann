/* compute_molecular_partition_probability.c
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
 * 
 *
 *  Created on: May 2, 2013
 *      Author:  Dennis G. Thomas
 *      Modified by Doug Baxter on May 31,2013.
 ********************************************************************************/

#include "boltzmann_structs.h"


#include "compute_molecular_partition_probability.h"

int compute_molecular_partition_probability(struct state_struct *state){
  /*
    Compute the partition_funcs and probabilities arrays (length unique_molecules)
    Called by: compute_standard_energies
    Calls:     exp,fopen,fprintf,fflush
  */

  struct molecule_struct *cur_molecules;
  char *molecules_text;
  char *molecule;

  double* molecule_dg0tfs;
  
  /*
    double  *partition_funcs; //len = unique_molecules
  */
  double  *probabilities; /* len = unique_molecules */
  
  double min_molecule_dg0tf;
  double scaling;
  double q;
  double temp_kelvin;
  double m_rt;
  double m_r_rt;

  int success;
  int i;

  int ci;
  int nu_molecules;

  int print_output;
  int padi;

  FILE *lfp;
  FILE *efp;

  success = 1;

  nu_molecules             = state->nunique_molecules;
  cur_molecules            = state->sorted_molecules;
  molecules_text           = state->molecules_text;

  min_molecule_dg0tf       = state->min_molecule_dg0tf;
  molecule_dg0tfs          = state->molecule_dg0tfs;
  /*
  partition_funcs          = state->molecular_partition_funcs;
  */
  probabilities            = state->molecule_probabilities;
  

  temp_kelvin              = state->temp_kelvin;
  m_rt                     = state->m_rt;
  m_r_rt                   = state->m_r_rt;

  print_output             = state->print_output;
  lfp                      = state->lfp;

  q       		   = 0.0;
  scaling 		   = min_molecule_dg0tf * m_r_rt;

  if (print_output){
    if (lfp) {
      fprintf(lfp, 
	      "Output from compute_molecular_partition_probability.c: \n\n");
      fprintf(lfp,"m_r_rt = %le, min_molecule_dg0tf = %le, scaling = %le\n",
	      m_r_rt,min_molecule_dg0tf,scaling);
      
    }
  }
  for (i=0;i<nu_molecules;i++) {
    if (print_output) {
      if (lfp) {
	fprintf(lfp,
		"deltaG0_tf of molecule at index %i = %f.\n",
		i,molecule_dg0tfs[i]);
      }
    }
    /*
      partition_funcs[i] = exp(-1.0*molecule_dg0tfs[i]/(IDEAL_GAS_CONST_KJ_PER_KELVIN_MOL*temp_kelvin));
      computing this sum may be a bit unstable.
    */
    probabilities[i] = exp((molecule_dg0tfs[i]*m_r_rt) - scaling);
    q = q + probabilities[i];
  }

  /*
    fes->molecular_partition_funcs_sum = q;
  */
  if (print_output) {
    if (lfp) {
      fprintf(lfp,"total partition function (q) = %f.\n",q);
    }
  }
  if (q <= 0.0) {
    fprintf(stderr,"compute_molecular_partition_probability: Error: "
	    " non-positive sum of partitions, q = %le\n",q);
    success = 0;
  }
  if (success) {
    q = 1.0/q;
    for (i=0;i<nu_molecules;i++) {
      molecule    = (char *)&molecules_text[cur_molecules->string];
      if (print_output) {
	if (lfp) {
	  fprintf(lfp,"partition function of molecule %s = %f.\n",
		  molecule,probabilities[i]);
	}
      }
      probabilities[i] = probabilities[i]*q;
      if (print_output) {
	if (lfp) {
	  fprintf(lfp,"probability of molecule %s = %f.\n",
		  molecule,probabilities[i]);
	  fprintf(lfp,"\n");
	}
      }
      cur_molecules += 1; /* Caution address arithmetic. */
    }
    if (print_output) {
      if (lfp) {
	fprintf(lfp,"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
      }
    }
  }
  return (success);
}
