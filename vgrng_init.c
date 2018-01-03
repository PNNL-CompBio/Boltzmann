/* vgrng_init.c
*******************************************************************************
Adapted from the MSpolygraph code vgrng_init.c:
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
#include "boltzmann_structs.h"
#include "vgrng.h"
#include "vgrng_init.h"

int64_t vgrng_init(struct vgrng_state_struct *vgrng_state, int n) {
  int64_t mask;
  int64_t sv;
  int i;
  /*
    Initialize the pseudo random number generator state, 
    pseudo random (uniformly distributed) 31 bit integers.
    Called by: boltzmann_init.
  */
  mask = (((int64_t)1)<<32) - ((int64_t)1);
  vgrng_state->mask = mask;
  vgrng_state->fib_constant     = ((int64_t)1597258441) & mask;
  vgrng_state->lcg_constant     = ((int64_t)2584418167) & mask;
  vgrng_state->fib_history[0]   = ((int64_t)6765109461) & mask;
  vgrng_state->fib_history[1]   = ((int64_t)1771128657) & mask;
  vgrng_state->lcg_history      = ((int64_t)4636875025) & mask;
  vgrng_state->fib_cur_ptr      = 1;
  sv = 0;
  for (i=0;i<n;i++) {
    sv += (vgrng(vgrng_state) & mask);
  }
  vgrng_state->uni_multiplier = 1.0/((double)mask);
  return(sv);
}

