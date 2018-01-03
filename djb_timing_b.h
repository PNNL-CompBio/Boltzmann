/* djb_timing_b.h
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
/*
  This file is used in building the .o files that go into
  the djb_timing.a library. Included by timing_init.c timing_print.c
  timing_start.c and timing_stop.c
*/
#ifndef _djb_timing_b_h
#define _djb_timing_b_h

#include <stdio.h>
#include <sys/types.h>

#ifdef	__cplusplus
extern "C"
{
#endif
  /*
    extern declarations for timing functions.
    may want to keep the  pargs and doprnt routines?
  */
#include <time.h> 
struct clock_data {
  int64_t acc_time;
  int64_t local_start_count;
  int64_t local_stop_count;
  int64_t call_count;
  char *tagcs;
  int  indent,dummy;
}
;
#define TIMING_MAX_ENTRIES 1024
#define TIMING_MAX_TAG_LEN   31
#define TIMINGI_BUF_LEN     120
#define TIMING_MAX_INDENT    10  

struct timing_struct {
  int64_t global_start_count;
  int64_t global_stop_count;
  struct clock_data timings[TIMING_MAX_ENTRIES];
  /*
  struct clock_data selftime[4];
  */
  double secs_per_clock;
  int64_t initialized;
  int64_t ticks_per_envelopes;
  int index;
  int start_stop_index;
  int init_index;
  int envelope_index;
  int user_index;
  int alpha_present;
  int dummy;
  char timings_alpha[TIMING_MAX_ENTRIES*TIMING_MAX_TAG_LEN+TIMING_MAX_ENTRIES];
}
;

#define TIMING_ITERATES 1000000
#define TIMING_NEWTIME -1
#define TIMING_INITIALIZED 42

#include "djb_timing.h"

#ifdef	__cplusplus
}
#endif
#endif

#define DJB_TIMING_START_STOP   0
#define DJB_TIMING_INIT         1
#define DJB_TIMING_ENVELOPE     2
#define DJB_TIMING_MIN          3
#define DJB_TIMING_MAX       1023
