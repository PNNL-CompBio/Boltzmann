/* alloc5.c
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

#include "alloc5.h"
int alloc5(int64_t num_cpds, 
	   int64_t pi_len, 
	   int64_t align_len, 
	   int64_t *usage,
	   struct formation_energy_struct **formation_energies_p) {
  /*
    Allocate space for the formation_energies struct, and its
    pseudoisomers, pseudoisomer_strings, sorted_pseudoisomer fields:
    space to hold data read from state->pseudoisomer_dg0_file.

    num_cpds is the number of pseudoisomers.

    pi_len is the length of the pseudoisomer_strings character array.

    align_len is the alignment value for strings in pseudoisomer_strings.

    *usage is incremented by the amount of space allocated;

    formation_energies_p is returned as a pointer to the allocated
    foromation_energy struct. 

    Called by: compute_standard_energies
    Calls:     calloc, fprintf (intrinsic)
  */
  struct formation_energy_struct fe_s;
  struct formation_energy_struct *formation_energies;
  struct pseudoisomer_struct pi_s;
  struct pseudoisomer_struct *pseudoisomers;
  char *pseudoisomer_strings;
  int64_t ask_for;
  int64_t one_l;
  int64_t enum_cpds;
  int64_t align_mask;
  int success;
  int padi;

  success = 1;
  one_l      = (int64_t)1;
  align_mask = align_len - one_l;
  /*
    Allocate space for the formation_energy struct.
  */
  ask_for = (int64_t)sizeof(fe_s);
  *usage += ask_for;
  formation_energies = (struct formation_energy_struct *)calloc(one_l,ask_for);
  
  if (formation_energies == NULL) {
    fprintf(stderr,"alloc5: Error, unable to allocate %lld bytes of space for "
	    "formation_energy struct\n",ask_for);
    fflush(stderr);
    success = 0;
  }
  if (success) {
    *formation_energies_p = formation_energies;
    /*
      Allocate space for the pseudoisomer dg0f data.
    */
  
    ask_for = num_cpds * ((int64_t)sizeof(pi_s));
    *usage   += ask_for;
    pseudoisomers = (struct pseudoisomer_struct *)calloc(one_l,ask_for);
    if (pseudoisomers == NULL) {
      fprintf(stderr,"alloc5: Error, unable to allocate %lld bytes of space "
	      "for pseudoisomer data \n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  if (success) {
    formation_energies->pseudoisomers = pseudoisomers;
    ask_for = pi_len + ((align_len - (pi_len & align_mask)) & align_mask);
    *usage += ask_for;
    pseudoisomer_strings = (char *)calloc(one_l,ask_for);
    if (pseudoisomer_strings == NULL) {
      fprintf(stderr,"alloc5: Error, unable to allocate %lld bytes of space "
	      "for pseudoisomer_strings data \n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  if (success) {
    formation_energies->pseudoisomer_strings = pseudoisomer_strings;
  }
  /*
  if (success) {
    enum_cpds = num_cpds + (num_cpds & 1);
    ask_for = enum_cpds * ((int64_t)sizeof(int));
    *usage += ask_for;
    sorted_pseudoisomers_indices = (int*)calloc(one_l,ask_for);
    if (sorted_pseudoisomer_indices == NULL) {
      fprintf(stderr,"alloc5: Error, unable to allocate %lld bytes of space "
	      "for sorted_pseudoisomer_indices.\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  if (success) {
    formation_energies->sorted_pseudoisomer_indices = 
      sorted_pseudoisomer_indices;
    ask_for = (enum_cpds + enum_cpds) * ((int64_t)sizeof(int));
    *usage += ask_for;
    pseudoisomer_scratch = (int*)calloc(one_l,ask_for);
    if (pseudoisomer_scratch == NULL) {
      fprintf(stderr,"alloc5: Error, unable to allocate %lld bytes of space "
	      "for pseudoisomer_scratch.\n",ask_for);
      fflush(stderr);
      success = 0;
    } else {
      formation_energies->pseudoisomer_scratch = pseudoisomer_scratch;
    }
  }
  */
  return(success);
}
