/* blank_to_dash.c
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
#include "blank_to_dash.h"

void blank_to_dash (int sl, char *in_string, char *out_string) {
  /*
    Convert blanks to dashes in a string.
    Called by: parse_pseudoisomer_dg0f_file
    Leaf routine.
  */
  /*
  int64_t add_con_in;
  int64_t add_con_out;
  int64_t i64_mask;
  int64_t *oi;
  int64_t *ii;
  int64_t blank_start;
  int64_t dash_delta;
  int64_t mk_start;
  int64_t dash;
  int64_t mk;
  int64_t cc;
  int64_t blank;
  int64_t adj;
  int64_t v;
  int lim;
  int j;
  */
  int ic;
  int istart;
  int i;
  int pad1;
  /*

  i64_mask = (int64_t)7;
  blank_start = ((int64_t)32) << 56;
  dash_delta  = ((int64_t)13) << 56;
  mk_start = i64_mask << 56;
  add_con_in = (int64_t)in_string;
  add_con_out = (int64_t)out_string;
  */
  /*
    If character strings are both eight byte aligned work
    on eight byte integers as much as possible.
  */
  /*
  if (((i64_mask & add_con_in) | (i64_mask & (add_con_out))) == 0) {
    ii = (int64_t *)in_string;
    oi = (int64_t *)out_string;
    istart = sl - (sl & i64_mask);
    for (i = 0; i<istart; i+= 8) {
      mk = mk_start;
      blank = blank_start;
      dash  = dash_delta;
      v  = *ii;
      for (j=0;j<8;j++) {
	cc  = v & mk;
	adj = 0 - (cc == blank);
	adj = adj & dash;
	v   += adj;
	mk  = mk >> 8;
	blank = blank >> 8;
	dash  = dash >> 8;
      }
      *oi = v;
      ii += 1;
      oi += 1;
    }
  } else {
    istart = 0;
  }
  */
  istart = 0;
  for (i=istart;i<sl;i++) {
    ic = (int)in_string[i];
    if (ic == 32) {
      ic = 45;
    }
    out_string[i] = (char)ic;
  }
}

