/* alloc3.c
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

#include "alloc3.h"

int alloc3(struct state_struct *state) {
  /*
    Allocate space for the following vector fields for reading in
    information about species concentrations, experimental values
    and user values.
    compartment_pointers   	     (nunique_compartments + 1) 
    current_counts         	     (nunique_molecules)
    bndry_flux_counts      	     (nunique_molecules)
    count_to_conc          	     (nunique_molecules)
    conc_to_count          	     (nunique_molecules)
    dg0s                   	     (number_reactions)  
    ke                     	     (number_reactions) 
    kss                    	     (number_reactions)
    kssr                   	     (number_reactions)
    kss_eval               	     (nunique_molecules)
    kss_uval               	     (nunique_molecules)
    dg0tfs                 	     (nunique_molecules)
    molecule_probabilities 	     (nunique_molecules)
    molecule_chemical_potentials     (nunique_molecules)
    Called by: boltzmann_init
    Calls:     calloc, fprintf, fflush (intrinsic)
  */
  struct molecules_matrix_struct sms;
  struct molecules_matrix_struct *molecules_matrix;
  int64_t ask_for;
  int64_t one_l;
  int64_t usage;
  int64_t align_len;
  int64_t align_mask;
  int nu_molecules;
  int max_molecule_len;

  int success;
  int nrxns;

  int nzr;
  int max_compartment_len;
  success = 1;
  one_l      		  = (int64_t)1;
  usage      		  = state->usage;
  align_mask 		  = state->align_mask;
  align_len  		  = state->align_len;
  nu_molecules            = (int)state->nunique_molecules;
  nzr                     = (int)state->number_molecules;
  nrxns                   = (int)state->number_reactions;
  max_molecule_len        = (int)state->max_molecule_len + 1;
  max_compartment_len     = (int)state->max_compartment_len + 1;
  /*
    Allocate space for compartment pointers in the sorted molecules list -
    length is unique_compartments + 1;
  */
  if (success) {
    ask_for = ((int64_t)state->nunique_compartments + one_l) * 
      ((int64_t)sizeof(int64_t));
    state->compartment_ptrs = (int64_t*)calloc(one_l,ask_for);
    if (state->compartment_ptrs == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "compartment_ptrs field\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }

  if (success) {
    /*
      Allocate space for the current counts buffer.
    */
    ask_for = ((int64_t)nu_molecules) * ((int64_t)sizeof(double));
    usage += ask_for;
    state->current_counts = (double *)calloc(one_l,ask_for);
    if (state->current_counts == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "current_counts field\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  if (success) {
    /*
      Allocate space for the boudary flux concentrations buffer.
    */
    ask_for = ((int64_t)nu_molecules) * ((int64_t)sizeof(double));
    usage += ask_for;
    state->bndry_flux_counts = (double *)calloc(one_l,ask_for);
    if (state->bndry_flux_counts == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "bndry_flux_counts field\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  if (success) {
    /*
      Allocate space for the count_to_conc vector.
    */
    ask_for = ((int64_t)nu_molecules) * ((int64_t)sizeof(double));
    usage += ask_for;
    state->count_to_conc = (double *)calloc(one_l,ask_for);
    if (state->count_to_conc == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "count_to_conc field\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }    
  if (success) {
    /*
      Allocate space for the conc_to_count vector.
    */
    ask_for = ((int64_t)nu_molecules) * ((int64_t)sizeof(double));
    usage += ask_for;
    state->conc_to_count = (double *)calloc(one_l,ask_for);
    if (state->conc_to_count == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "conc_to_count field\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }    

  /*
    Allocate space for the change in Gibb's free energy.
  */
  if (success) {
    ask_for = ((int64_t)nrxns) * ((int64_t)sizeof(double));
    usage += ask_for;
    state->dg0s = (double *)calloc(one_l,ask_for);
    if (state->dg0s == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "state->dg0s field.\n",ask_for);
      fflush(stderr);
      success = 0;
    } 
  }
  /*
    Allocate space for the reaction equilibrium coefficients.
  */
  if (success) {
    ask_for = ((int64_t)nrxns) * ((int64_t)sizeof(double));
    usage += ask_for;
    state->ke = (double *)calloc(one_l,ask_for);
    if (state->ke == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "state->ke field.\n",ask_for);
      fflush(stderr);
      success = 0;
    } 
  }
  /*
    Allocate space for the steady state reaction equilibrium coefficients
    adjustments..
  */
  if (success) {
    ask_for = ((int64_t)nrxns) * ((int64_t)sizeof(double));
    /*
      We also need to store a kss for each reverse reaction as
      they are not reciprocals.
    */
    ask_for += ask_for;
    usage += ask_for;
    state->kss = (double *)calloc(one_l,ask_for);
    if (state->kss == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "state->kss field.\n",ask_for);
      fflush(stderr);
      success = 0;
    } else {
      state->kssr = (double *)&(state->kss[nrxns]);
    }
  }
  /*
    Allocate space for experimental concentration.
  */
  if (success) {
    ask_for = ((int64_t)nu_molecules) * ((int64_t)sizeof(double));
    usage += ask_for;
    state->kss_e_val = (double *)calloc(one_l,ask_for);
    if (state->kss_e_val == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "state->kss_e_val field.\n",ask_for);
      fflush(stderr);
      success = 0;
    } 
  }
  /*
    Allocate space for user concentrataion
  */
  if (success) {
    ask_for = ((int64_t)nu_molecules) * ((int64_t)sizeof(double));
    usage += ask_for;
    state->kss_u_val = (double *)calloc(one_l,ask_for);
    if (state->kss_u_val == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "state->kss_u_val field.\n",ask_for);
      fflush(stderr);
      success = 0;
    } 
  }
  /*
    Allocate space for molecule_dg0tfs
  */
  if (success) {
    ask_for = ((int64_t)nu_molecules) * ((int64_t)sizeof(double));
    usage += ask_for;
    state->molecule_dg0tfs = (double *)calloc(one_l,ask_for);
    if (state->molecule_dg0tfs == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "state->molecule_dg0tfs field.\n",ask_for);
      fflush(stderr);
      success = 0;
    } 
  }
  /*
    Allocate space for molecule_probabilities
  */
  if (success) {
    ask_for = ((int64_t)nu_molecules) * ((int64_t)sizeof(double));
    usage += ask_for;
    state->molecule_probabilities = (double *)calloc(one_l,ask_for);
    if (state->molecule_probabilities == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "state->molecule_probabilities field.\n",ask_for);
      fflush(stderr);
      success = 0;
    } 
  }
  /*
    Allocate space for molecule_chemical_potentials
  */
  if (success) {
    ask_for = ((int64_t)nu_molecules) * ((int64_t)sizeof(double));
    usage += ask_for;
    state->molecule_chemical_potentials = (double *)calloc(one_l,ask_for);
    if (state->molecule_chemical_potentials == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "state->molecule_chemical_potentials field.\n",ask_for);
      fflush(stderr);
      success = 0;
    } 
  }
  state->usage = usage;
  return(success);
}
