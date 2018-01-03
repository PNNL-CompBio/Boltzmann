/* vgrng.c
*******************************************************************************
MSPolygraph

Pacific Northwest National Laboratory, Richland, WA 99352.

MSPolygraph is a mass spectrogram analysis tool that identifies likely 
matching candidate peptides from a reference database.

Copyright (c) 2010 Battelle Memorial Institute.

Publications based on work performed using the software should include 
the following citation as a reference:

    Cannon WR, KH Jarman, BM Webb-Robertson, DJ Baxter, CS Oehmen, 
    KD Jarman, A Heredia-Langner, GA Anderson, and KJ Auberry.
    A Comparison of Probability and Likelihood Models for Peptide 
    Identification from Tandem Mass Spectrometry Data, 
    J. Proteome Res., 2005 4(5):1687-1698

Licensed under the Educational Community License, Version 2.0 (the "License"); 
you may not use this file except in compliance with the License. 
The terms and conditions of the License may be found in 
ECL-2.0_LICENSE_TERMS.TXT in the directory containing this file.
        
Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR 
CONDITIONS OF ANY KIND, either express or implied. See the License for the 
specific language governing permissions and limitations under the License.
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <float.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <errno.h>


#include "vgrng_state_struct.h"

#include "vgrng.h"

int64_t vgrng(struct vgrng_state_struct *vgrng_state) {
  int64_t a,b,c,d,e,mask;
  int64_t threec,twelvec,thirteenc,f;
  int i;
  i = vgrng_state->fib_cur_ptr;
  a = vgrng_state->fib_history[i];
  b = vgrng_state->fib_history[1-i];
  c = vgrng_state->lcg_history;
  d = vgrng_state->fib_constant;
  e = vgrng_state->lcg_constant;
  mask = vgrng_state->mask;
  b = (a+b+d) & mask;
  i = 1-i;
  vgrng_state->fib_history[i] = b;
  vgrng_state->fib_cur_ptr = i;
  /*
  c = (c * vgrng_state->lcg_multplier + e) & mask;
  */
  /*
    Hardwiring in the multiplier 69069 for now to
    remove the integer multiply though 
    I do not understand why chip manufacturers can't build a
    fast integer multiply using Napier's Bones
  */
  f = c << 16;
  threec = (c << 1) + c;
  twelvec = threec << 2;
  thirteenc = twelvec+c;
  f += (thirteenc << 8);
  f += (twelvec << 4);
  f += thirteenc;
  c = (f + e) & mask;
  vgrng_state->lcg_history = c;
  f = b^c;
  vgrng_state->xor_lcg_fib = f;
  return(f);
}
