/* compute_chemical_potential.c
******************************************************************************
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
 * compute_chemical_potential.c
 *
 *  Created on: May 2, 2013
 *      Author:  Dennis G. Thomas
 *  Adapted by doug baxter.5/30/2013
 *  Compute the chemical_potentials from the current concentrations and
 *  the probabilities vector.
 *
 *  Called by: compute_standard_energies
 *  Calls:     log
 
*/

#include "boltzmann_structs.h"
#include "compute_chemical_potential.h"

int compute_chemical_potential(struct formation_energy_struct *fes){

  double  *probabilities; /*len = unique_molecules*/
  double  *chemical_potentials; /*len = unique_molecules */
  double *current_counts;	/*len = unique_molecules */
  double sum;
  double sump1;
  double m_rt;

  int success;
  int i;
  int nu_molecules;
  int pad_i;

  success = 1;
  

  nu_molecules        = fes->nunique_molecules;
  probabilities       = fes->molecule_probabilities;
  chemical_potentials = fes->molecule_chemical_potentials;
  current_counts      = fes->current_counts;

  /*
    -RT is computed in read_params and stored in state->m_rt and 
    in turn in fes->m_rt;
  */

  m_rt = fes->m_rt;

  sum = 0.0;
  for (i=0;i<nu_molecules;i++) {
    sum = sum + current_counts[i];
  }
  sump1 = sum + 1.0;
  for (i=0;i<nu_molecules;i++) {
    chemical_potentials[i] = m_rt *
      log((sump1*probabilities[i])/(current_counts[i]+1.0));
  }
  return (success);
}
