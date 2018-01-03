/* upcase.c
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
#include "upcase.h"

void upcase (int sl, char *in_string, char *out_string) {
  /*
    Convert lower case to upper case in a string that may
    have non-alphabetic characters included.
    Called by: parse_reactions_file, parse_pseudoisomer_dg0f_file
    Leaf routine.
  */
  /*
  int64_t add_con_in;
  int64_t add_con_out;
  int64_t i64_mask;
  int64_t *oi;
  int64_t *ii;
  int64_t hi_start;
  int64_t lo_start;
  int64_t uc_start;
  int64_t mk_start;
  int64_t hi;
  int64_t lo;
  int64_t uc;
  int64_t mk;
  int64_t cc;
  int64_t adj;
  int64_t v;
  int lim;
  int j;
  */
  int i;
  int ic;
  int istart;
  int pad1;
  /*
  i64_mask = (int64_t)7;
  hi_start = ((int64_t)123) << 56;
  lo_start = ((int64_t)96) << 56;
  uc_start = ((int64_t)32) << 56;
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
      uc = uc_start;
      hi = hi_start;
      lo = lo_start;
      v  = *ii;
      for (j=0;j<8;j++) {
	cc  = v & mk;
	adj = 0 - ((cc > lo) & (cc < hi));
	adj = adj & uc;
	v   -= adj;
	mk  = mk >> 8;
	uc  = uc >> 8;
	hi  = hi >> 8;
	lo  = lo >> 8;
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
    if ((ic > 96) && (ic < 123)) {
      ic = ic - 32;
    }
    out_string[i] = (char)ic;
  }
}

