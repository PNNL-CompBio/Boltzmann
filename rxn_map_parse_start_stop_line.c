/* rxn_map_parse_start_stop_line.c
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

#include "count_nws.h"
#include "count_ws.h"
#include "find_colon.h"
#include "compartment_lookup.h"
#include "molecules_lookup.h"

#include "rxn_map_parse_start_stop_line.h"
int rxn_map_parse_start_stop_line(struct state_struct *state, 
				  int *first,
				  int *last) {
  /*
    Parse a line of a start_stop file setting *first to the index
    of the first molecule, and *last to the index of the second molecule.
    the start_stop_fp is taken from state->conc_fp and assumed to
    have been opened in rxn_map_init.
    Returns 1 on success, 0 on failure.
    Called by: rxn_map
    Calls:     count_nws, 
               count_ws, 
               find_colon, 
	       compartment_lookoup,
               molecules_lookup,
	       fgets, fprintf, fflush
  */
  char *line_buffer;
  char *lbp;
  char *start_mol;
  char *start_cmpt;
  char *stop_mol;
  char *stop_cmpt;

  int  buff_size;
  int  pos;

  int  sci;
  int  eci;

  int  start_len;
  int  stop_len;

  int  success;
  int  padi;

  FILE *start_stop_fp;
  start_stop_fp = state->conc_fp;
  line_buffer   = state->param_buffer;
  success = 1;
  *first = -1;
  *last  = -1;
  buff_size = (int)state->max_param_line_len;
  if (start_stop_fp == NULL) {
    success = 0;
  } 
  if (success) {
    lbp = fgets(line_buffer,buff_size,start_stop_fp);
    if (lbp == NULL) {
      success = 0;
    }
  }
  if (success) {
    pos = count_ws(line_buffer);
    line_buffer += pos; /* Caution address arithmetic. */
    start_mol = line_buffer;
    start_len   = count_nws(start_mol);
    if (start_len == 0) {
      success = 0;
    }
  }
  if (success) {
    line_buffer += start_len; /* Caution address arithmetic. */
    pos = count_ws(line_buffer);
    line_buffer += pos; /* Caution address arithmetic. */
    if (pos == 0) {
      success = 0;
    }
  }
  if (success) {
    stop_mol = line_buffer;
    stop_len = count_nws(stop_mol);
    if (stop_len == 0) {
      success = 0;
    }
  }
  if (success) {
    /*
      Terminate start and stop molecule strings.
    */
    start_mol[start_len] = '\0';
    stop_mol[stop_len]   = '\0';
    /*
      Look for a compartment for start molecule.
    */
    pos = find_colon(start_mol);
    if (pos < 0) {
      sci = 0;
    } else {
      start_cmpt = (char *)&start_mol[pos+1];
      start_mol[pos] = '\0';
      sci = compartment_lookup(start_cmpt,state);
    }
    *first = molecules_lookup(start_mol,sci,state);
    if (*first < 0) {
      if (sci <= 0) {
	fprintf(stderr,"rxn_map_parse_start_stop_line: Error, molecule %s not found\n",start_mol);
      } else {
	fprintf(stderr,"rxn_map_parse_start_stop_line: Error, molecule %s:%s not found\n",start_mol,start_cmpt);
      }
      fflush(stderr);
      success = 0;
    }
  }
  if (success) {
    pos = find_colon(stop_mol);
    if (pos < 0) {
      eci = 0;
    } else {
      stop_cmpt = (char*)&stop_mol[pos+1];
      stop_mol[pos] = '\0';
      eci = compartment_lookup(stop_cmpt,state);
    }
    *last = molecules_lookup(stop_mol,eci,state);
    if (*last < 0) {
      if (eci <= 0) {
	fprintf(stderr,"rxn_map_parse_start_stop_line: Error, molecule %s not found\n",stop_mol);
      } else {
	fprintf(stderr,"rxn_map_parse_start_stop_line: Error, molecule %s:%s not found\n",stop_mol,stop_cmpt);
      }
      fflush(stderr);
      success = 0;
    }
  }
  return (success);
}
    
