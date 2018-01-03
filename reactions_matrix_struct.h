/* reactions_matrix_struct.h 
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
#ifndef _REACTIONS_MATRIX_STRUCT_H_ 
#define _REACTIONS_MATRIX_STRUCT_H_  1
struct reactions_matrix_struct {
  /*
    This structure is laid out in the form a row-compressed
    sparse matrix structure, but with 4 sets of "column numbers"
    One could also view it as each entry has three components,
    a molecule index, a compartment index, a coefficient,  and a text pointer.
    The rxn_prtrs vector is of length (number_reactions + 1).
    rxn_ptrs[i] points to the first "element" position for reaction
    number i (i=0;number_reactions - 1) and rxn_ptrs[number_reactions] is
    1 more than the total number of molecules in all the reactions
    (without regard to uniqueness). The the data for reaction i, is in
    the elements of the four vectors molecules_indices, compartment_indices,
    coefficients, and text from element rxn_ptr[i] through element 
    rxn_ptr[i+1]-1.

    molecules_indices, coefficients and text fields are first filled in
    by parse_side_line called by parse_reactions_file.

    unique_compartments fills in the compartment_indices fields.
    These fields are all allocated in alloc2.

  */

  int64_t *rxn_ptrs;
  /*
    molecules_indices indicate which molecules are present as reactants
    and products.
  */
  int64_t *molecules_indices;
  /*
    Compartment indices indicate which compartments the molecules
    are in.
  */
  int64_t *compartment_indices;
  /*
    coefficients are negative for reactants, positive for products.
  */
  int64_t *coefficients;
  /*
    Offsets into the molecules_text for the molecules name,
    used to fill in the molecules indices after molecules have been sorted.
  */
  int64_t *text;
  /*
    Since the coefficients for the solvent molecule are 
    zeroed out to run the simulation we need to be
    able to recover them. Sinced we only allow one solvemt molecule
    there will be at most num_rxns solvent coefficients.
  */
  int64_t *solvent_coefficients;
}
;
#endif
