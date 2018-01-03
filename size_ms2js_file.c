/* size_ms2js_file.c
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

#include "size_ms2js_file.h"
int size_ms2js_file(struct state_struct *state,
		    int64_t *num_modelseed_ids,
		    int64_t *length_ms2js_strings) {
  /*
    Determine the number of lines in the ms2js_file, and the
    total size of the file so as to set the num_modelseed_ids argument,
    and the length_ms2js_strings argument.
    This routine uses a system call to wc to get this information instead
    of actually reading the file. For more detailed size estimate we
    could mimic the flavor of size_rxns_file.c
    For each entry there are 2 string fields, 
    for padding we need 2*align_len * num_modelseed_ids more characters in the
    ms2js_strings field.


    Called by: boltzmann_init.
    Calls:     fpintf,fflush,fopen,fgets,sscanf,fclose,system
  */
  char *command;
  char *ms2js_file;
  char buffer[1048];
  int64_t lc;
  int64_t wc;
  int64_t cc;
  int success;
  int name_len;

  int buff_len;
  int ns;

  FILE *ms2js_size_fp;
  FILE *lfp;
  success = 1;
  command = (char*)&buffer[0];
  ms2js_file = state->ms2js_file;
  lfp        = state->lfp;
  name_len = (int)strlen(ms2js_file);
  buff_len = 1047;
  if (name_len + 23 > 1048) { 
    if (lfp) {
      fprintf(lfp,"size_ms2js_file: Error ms2js_file name is too long.\n");
      fflush(lfp);
    }
    success = 0;
  } else {
    sprintf(command,"wc %s > _ms2js_wc_output_",ms2js_file);
    system(command);
    ms2js_size_fp = fopen("_ms2js_wc_output_","r");
    if (ms2js_size_fp == NULL) {
      if (lfp) {
	fprintf(lfp,"size_ms2js_file: Error, unable to open _ms2js_wc_output_ "
	      "for sizing ms2js_file\n");
	fflush(lfp);
      }
      success = 0;
    } else {
      fgets(command,buff_len,ms2js_size_fp);
      ns = sscanf(command,"%lld %lld %lld",&lc,&wc,&cc);
      if (ns !=3) {
	if (lfp) {
	  fprintf(lfp,
		  "size_ms2js_file: Error, reading output of wc command\n");
	  fflush(lfp);
	}
	success = 0;
      } else {
	*num_modelseed_ids = lc;
	*length_ms2js_strings = cc  + 
	  ((state->align_len << 1) * lc);
	fclose(ms2js_size_fp);
	sprintf(command,"/bin/rm -f _ms2js_wc_output_");
	system(command);
      }
    }
  }
  return(success);
}
