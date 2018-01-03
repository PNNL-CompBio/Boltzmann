/* read_compartment_sizes.c
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

#include "read_compartment_sizes.h"

int read_compartment_sizes(struct state_struct *state) {
  /*
    Read in the sizes of compartments.
    Called by: species_init
    Calls:     fopen, fprintf, fflush, fgets, fclose
  */
  struct compartment_struct *compartments;
  struct compartment_struct *compartment;
  char *compartment_file;
  char *cmpts_buffer;
  char *compartment_name;
  char *line;
  char *fgp;
  char *vol_string;
  int64_t cmpts_buff_len;
  double volume;
  double units_avo;
  double multiplier;
  int success;
  int padi;

  int ci;
  int nb;

  int nc;
  int ns;

  FILE *cmpt_fp;
  FILE *lfp;
  FILE *error_fp;
  success = 1;
  compartment_file = state->compartment_file;
  cmpts_buff_len   = (int)state->max_param_line_len;
  cmpts_buffer     = state->param_buffer;
  lfp              = state->lfp;
  compartments     = state->sorted_cmpts;
  units_avo        = state->conc_units * state->avogadro;
  if (lfp) {
    error_fp = lfp;
  } else {
    error_fp = stderr;
  }
  cmpt_fp = fopen(compartment_file,"r");
  if (cmpt_fp == NULL) {
    fprintf(error_fp,"read_compartment_sizes: Error, unable to open file %d\n",
	    compartment_file);
    fflush(error_fp);
    success = 0;
  }
  if (success) {
    while (!feof(cmpt_fp)) {
      fgp = fgets(cmpts_buffer,cmpts_buff_len,cmpt_fp);
      compartment_name = cmpts_buffer;
      /*
	Skip over leading white space .
      */
      nb = count_ws(cmpts_buffer);
      compartment_name = cmpts_buffer + nb; /* Caution address arithmetic */
      nc = count_nws(compartment_name);
      compartment_name[nc] = '\0';
      vol_string = (char*)(&compartment_name[nb+1]);
      /*
	remove leading white space.
      */
      nb = count_ws(vol_string);
      vol_string += nb; /* Caution address arithmetic */
      ns = sscanf(vol_string,"%le",&volume);
      if ((ns < 1) || (volume <= 0.0)) {
	volume = state->default_volume;
	fprintf(error_fp,
		"read_compartment_sizes: Warning invalid volume for "
		"compartment %s, using %le\n",
		compartment_name,volume);
	fflush(error_fp);
      }
      ci = compartment_lookup(compartment_name,state);
      if (ci < 0) {
	fprintf(error_fp,"read_compartment_sizes: Warning, compartment %s "
		"not listed in reaction file.\n");
	fflush(error_fp);
      } else {
	compartment = (struct compartment_struct *)&compartments[ci];
	compartment->volume = volume;
	compartment->recip_volume = 1.0/volume;
	multiplier = units_avo * volume;
	compartment->conc_to_count = multiplier;
	compartment->count_to_conc = 1.0/multiplier;
      }
    }
    fclose(cmpt_fp);
  }
  return(success);
}
