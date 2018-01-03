/* sbml2bo_struct.h 
  
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
#ifndef __SBML2BO_STRUCT__
#define __SBML2BO_STRUCT__ 1

struct sbml2bo_struct {
  struct compartment_struct *sorted_compartments;
  struct compartment_struct *unsorted_compartments;
  char *sbml_file;
  char *concs_in_file;
  char *rxns_dat_file;
  char *cmpts_dat_file;
  char *log_file;
  char *compartment_text;
  double avogadro;
  double recip_avogadro;
  double default_comp_size;
  
  int  align_mask;
  int  alignment;
  int  file_name_len;
  int  num_files;

  int  num_reactions;
  int  num_species;

  int  num_cmpts;
  int  max_compartment_len;

  FILE *sbml_fp;
  FILE *concs_fp;
  FILE *rxns_fp;
  FILE *cmpts_fp;
  FILE *log_fp;
}
;
#endif
