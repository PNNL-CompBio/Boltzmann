/* pseudoisomer_struct.h 
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
#ifndef _PSEUDOISOMER_STRUCT_H_
#define _PSEUDOISOMER_STRUCT_H_ 1
struct pseudoisomer_struct {
  double z;  /* ionic strength */
  double dg0_f; /* free energy of formation*/
  /* offset of the json_cpd_name in the pseudoisomer_strings field */
  int64_t json_cpd_name; 
  /* offset of the kegg_id in the pseudoisomer_strings field */
  int64_t kegg_id;
  /* offset of the source string in the pseudoisomer_strings field */
  int64_t source;
  /* offset of the reference string in the pseudoisomer_strings field */
  int64_t ref;
  int  priority_level;
  int  nh; /* number of hydrogens */
}
;
#endif
