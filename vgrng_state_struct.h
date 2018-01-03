/* vgrng_state_struct.h
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
#ifndef _VGRNG_STATE_STRUCT_DEF_ 
#define _VGRNG_STATE_STRUCT_DEF_  1

struct vgrng_state_struct {
  /* 
     fib_history[1-cur_ptr] = (fib_history[cur_prtr] + fibhistory[1-curptr] + 
                               fib_constant) & mask; cur_ptr = 1-cur_ptr;
     lcg_history = (lcg_history*69069 + lcg_constant) & mask
  */
  int64_t fib_history[2]; /* history for Fibonacci generator*/
  int64_t lcg_history;    /* lcg history */
  int64_t fib_constant;   
  int64_t lcg_constant;
  int64_t lcg_multiplier;
  int64_t mask;
  int64_t fib_seed[2];
  int64_t lcg_seed;       
  int64_t xor_lcg_fib;
  double  uni_multiplier;
  int fib_cur_ptr;
  int padi;
  int64_t padl;
}
;
#endif
