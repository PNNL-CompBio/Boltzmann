/* istring_elem_struct.h 
  
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
*/
#ifndef __ISTRING_ELEM_STRUCT__
#define __ISTRING_ELEM_STRUCT__
struct istring_elem_struct {
  /*
    An istring_elem_struct has 4 fields. The first three are initially set by 
    the parse_side_line routine called by the parse_reactions_file routine.
    
        string  - pointer to a null terminated string that
	          contains the molecule or compartment string.

	m_index - If string points to a molecule, m_index
         	  is the ordinal position of the molecule in the
		  LEFT and RIGHT lines of the reactions.dat input file.
		  If the string is a compartment name m_index will be -1.

        c_index - If the string points to a compartment, c_index is the
	          ordinal position of the compartment name in the
		  COMPARTMENT, LEFT_COMPARTMENT, RIGHT_COMPARTMENT line, or
		  colon preceded fields of the LEFT and RIGHT lines in 
		  the reactions.dat input file.
		  If the string points to a molecule and the molecule
		  is in a compartment, c_index is set to the ordinal number
		  of the containing compartment, otherwise if the 
		  molecule is not in a compartment its c_index is -1.

	variable - An indicator set by the read_initial_concentrations
	          routine indicating whether or not this molecule is 
		  held fixed in concentration.
		  
	g_index	  The global index of the molecule or compartment, as set 
	          by the global reader.
	
  */
  int64_t string;
  int  m_index;
  int  c_index;
  int  variable;
  int  g_index;
}
;
#endif
