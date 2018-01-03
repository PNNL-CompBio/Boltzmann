/* ms2js_struct.h
*******************************************************************************
boltzmann
Author: Doug Baxter.

Pacific Northwest National Laboratory, Richland, WA 99352.

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
#ifndef _MS2JS_STRUCT_DEF_ 
#define _MS2JS_STRUCT_DEF_  1
struct ms2js_struct {
  int64_t num_modelseed_ids;
  int64_t length_ms2js_strings; 
  char    *ms2js_strings;
  char    **ms_ids;
  char    **js_ids;
}
;
#endif
