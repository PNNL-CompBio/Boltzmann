/* sbml_find_section.c
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

#include "count_ws.h"

#include "sbml_find_section.h"
int sbml_find_section(FILE *sbml_fp,char *sbml_buffer,int sbml_buffer_len,
		      char *section_name){
  /*
    Advance the file pointer sbml_fp to the first line after the section_name
    string exists. Returns a 1 if the section is found, 0 if not.
    modifies the file pointer for sbml_fp.
    Called by: parse_sbml, sbml_count_cmpts, sbml_count_species
    Calls:     feof,fgets,strlen,strncmp,count_ws
  */
  char *line;
  int in_section;
  int name_len;
  int success;
  int nb;
  success = 1;
  in_section = 0;
  name_len = strlen(section_name);
  
  while (!feof(sbml_fp) && (in_section == 0)) {
    fgets(sbml_buffer,sbml_buffer_len,sbml_fp);
    line = sbml_buffer;
    nb = count_ws(sbml_buffer);
    line += nb; /* Caution address arithmetic */
    if (strncmp(line,section_name,name_len) == 0) {
      in_section = 1;
    }
  }
  return (in_section);
}
