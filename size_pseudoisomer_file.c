/* size_pseudoisomer_file.c
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

#include "size_pseudoisomer_file.h"
int size_pseudoisomer_file(struct state_struct *state,
			   int64_t *num_pseudoisomers,
			   int64_t *length_pseudoisomer_strings) {
  /*
    Determine the number of entries in the pseudoisomer_file, and the
    total size of the file so as to set the num_pseudoisomers field,
    and the length_pseudoisomers_strings.
    This routine uses a system call to wc to get this information instead
    of actually reading the file. For more detailed size estimate we
    could mimic the flavor of size_rxns_file.c
    For each entry there are 4 string fields, and 4 numeric fields thus 
    for padding we need 4*align_len * num_pseudoisomers more characters in the
    pseudoisomers_strings field.
    This routine also sets the sue_pseudoisomers field of state to 1 if it
    could do a word count (wc) command on the pseudoisomer_file and 0 otherwise.

    Called by: compute_reaction_dg0fs_from_molecule_formation_energies.
    Calls:     fpintf,fflush,fopen,fgets,sscanf,fclose,system
  */
  char *command;
  char *dg0f_file;
  char buffer[1048];
  int64_t lc;
  int64_t wc;
  int64_t cc;
  int success;
  int command_len;
  int name_len;
  int buff_len;
  int ns;
  int np;
  FILE *dg0f_size_fp;
  success = 1;
  command = (char*)&buffer[0];
  dg0f_file = state->pseudoisomer_file;
  name_len = (int)strlen(dg0f_file);
  buff_len = 1047;
  if (name_len + 23 > 1048) { 
    fprintf(stderr,"size_pseudoisomer_file: Error pseudoisomer_file name is too long.\n");
    fflush(stderr);
    success = 0;
  } else {
    sprintf(command,"wc %s > _dg0f_wc_output_",dg0f_file);
    system(command);
    dg0f_size_fp = fopen("_dg0f_wc_output_","r");
    if (dg0f_size_fp == NULL) {
      fprintf(stderr,"size_pseudoiosmoer_file: Error, unable to open _dg0f_wc_output_ "
	      "for sizing pseudoisomer_dg0f_file\n");
      fflush(stderr);
      success = 0;
    } else {
      fgets(command,buff_len,dg0f_size_fp);
      ns = sscanf(command,"%ld %ld %ld",&lc,&wc,&cc);
      if (ns !=3) {
	fprintf(stderr,
		"size_pseudoisomer_file: Error, reading output of wc command\n");
	fflush(stderr);
	success = 0;
      } else {
	np = (int)lc - 1;
	*num_pseudoisomers = (int64_t)np;
	*length_pseudoisomer_strings = cc  + 
	  ((state->align_len << 2) * np);
	fclose(dg0f_size_fp);
	sprintf(command,"/bin/rm -f _dg0f_wc_output_");
	system(command);
      }
    }
  }
  return(success);
}
