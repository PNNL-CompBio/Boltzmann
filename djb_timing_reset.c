/* djb_timing_reset.c
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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <float.h>
#include "djb_timing_b.h"
/*
extern void djb_timing_start(int *index);
			 

extern void djb_timing_stop(int *index);
*/
extern struct timing_struct timing_data;

void djb_timing_reset() {
  /*
    Initialize the timings data structure.
  */
  int  i;
  for (i=0;i<TIMING_MAX_ENTRIES;i++) {
    (timing_data.timings)[i].acc_time=(int64_t)0;
    (timing_data.timings)[i].local_start_count=(int64_t)0;
    (timing_data.timings)[i].local_stop_count=(int64_t)0;
  }
}






