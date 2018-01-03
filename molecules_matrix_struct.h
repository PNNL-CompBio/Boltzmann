/* molecules_matrix_struct.h 
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
#ifndef __MOLECULES_MATRIX_STRUCT__ 
#define __MOLECULES_MATRIX_STRUCT__  1
struct molecules_matrix_struct {
  int64_t *molecules_ptrs;
  /*
    rxn_indices indicate which reactions the molecules are present in 
    as reactants or and products.
  */
  int64_t *reaction_indices;
  /*
    coefficients are negative for reactants, positive for products.
  */
  int64_t *coefficients;
}
;
#endif
