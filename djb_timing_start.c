/* djb_timing_start.c
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
#include <sys/types.h>
#include "djb_timing_b.h"

extern struct timing_struct timing_data;
void djb_timing_start(int index, int indent) {
  /*
  clock_t start;
  */
  int64_t start;
  int64_t mumble1;
  int64_t mumble2;
  int mumble_cnt;
  int i;
  int j;
  /* 
    Sanity check on input arguments.
  if (&timing_data == NULL) {
    fprintf(stderr,
	    " Warning: From timing_start: No timing_data allocated. \n");
    fprintf(stderr,
	    "          Please include \"timing.h\" in your main program.\n");
  } else {
    if (timing_data.initialized != TIMING_INITIALIZED) {
      fprintf(stderr,
	      "Warning: from timing_start: timing_data has not been initiallized. \n");
      fprintf(stderr,
	      "Please insert a call to TIMING_INIT(\"timingi.h\");\n");
      fprintf(stderr,"before the first call to TIMING_START in your code.\n");
    } else {
      if ((index < 0) || (index >= TIMING_MAX_ENTRIES)) {
	fprintf(stderr,
	       "Warning: from timing_start: invalid index ignored.\n");
	fprintf(stderr,
	       "         Indices must be >= 0 and < %d .\n",TIMING_MAX_ENTRIES);

	fprintf(stderr,
	       "         Index was %d .\n",index);
      } else {
  */
  start = itc_clock();
  i = index & DJB_TIMING_MAX;
  /*
    We will not count this call as inside itself.
  */
  timing_data.global_start_count += 1;
  ((timing_data.timings)[i]).local_start_count -= timing_data.global_start_count;
  ((timing_data.timings)[i]).local_stop_count -= timing_data.global_stop_count;
  ((timing_data.timings)[i]).call_count += 1;
  /* 
    do the mumble.
  */
  mumble_cnt = 7;
  mumble1 = timing_data.global_start_count;
  mumble2 = 69069L + i;
  for (j=0;j<mumble_cnt;j++) {
    mumble2 = mumble1 + mumble2;
    mumble1 = (mumble2 << 23) | (mumble2 >> 41);
  }
  if (mumble1 == 42) {
    fprintf(stderr," mumble \n");
    fflush(stderr);
  }
  ((timing_data.timings)[i]).acc_time -= start ;
  /*
      }
    }
  }       
  */    
}


