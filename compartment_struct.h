/* compartment_struct.h 
  
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
#ifndef __COMPARTMENT_STRUCT__
#define __COMPARTMENT_STRUCT__
struct compartment_struct {
  /*
    An molecule_struct has 8 fields. The first three are initially set by 
    the parse_side_line routine called by the parse_reactions_file routine.
    
        volume  - volume of the compartment, must be > 0
	recip_volume 1/volume

        string  - character offest into the molecules_text array 
	          to a null terminated string that contains the molecule 
		  or compartment string.
	          

        c_index - c_index is the
	          ordinal position of the compartment name in the
		  COMPARTMENT, LEFT_COMPARTMENT, RIGHT_COMPARTMENT line, or
		  colon preceded fields of the LEFT and RIGHT lines in 
		  the reactions.dat input file.

	g_index	  The global index of the compartment, as set 
	          by the global reader.
	
     For compartments only the string, c_index and g_index fields are used.

  */
  double volume;
  double recip_volume;
  double ntotal_exp;
  double ntotal_opt;
  /* 
     min_conc is really just count_to_conc.
  double min_conc;
  */
  double conc_to_count;
  double count_to_conc;

  int64_t string;

  int  c_index;
  int  g_index;
}
;
#endif
