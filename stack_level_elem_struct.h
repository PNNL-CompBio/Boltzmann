/* stack_level_elem_struct.h 
*******************************************************************************
Boltzmann

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
#ifndef __STACK_LEVEL_ELEM_STRUCT__
#define __STACK_LEVEL_ELEM_STRUCT__ 1
struct stack_level_elem_struct {
  int64_t first_stack_pos;
  int64_t cur_stack_pos;
  int64_t last_stack_pos;
}
;
#endif

