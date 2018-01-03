/* djb_timing_init.c
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

void djb_timing_init(char *timings_list_indx) {
  /*
    Initialize the timings data structure.
  */
  double  d_one;
  /*
  clock_t start; 
  */
  int64_t start;
  int64_t stopt;
  int64_t mumble1;
  int64_t mumble2;
  FILE *TIMINGI_FP;  
  int  noents,val,indent;
  int  i, j;
  int  scr_indent;
  int  mumble_cnt;
  char def_str[8];
  char desc[80];
  char timingi_buf[TIMINGI_BUF_LEN];
  start = itc_clock();
  d_one = 1.0;
  mumble_cnt = 10;
  /*
    Need to check that timing_data is a valid pointer.
  */
  /*
  fprintf (stdout, "djb_timing_init: &timing_data = %px\n",&timing_data);
  fflush(stdout);
  */
  /*
  if (timing_data == NULL) {
    fprintf(stderr," Error: in timing_init: timing_data pointer is NULL.\n");
    fprintf(stderr,
	    " You need to include \"timing.h\" in your main program.\n");
  } else {
  */
    if (timing_data.initialized != TIMING_INITIALIZED) {
      timing_data.initialized = TIMING_INITIALIZED;
      timing_data.global_start_count=0;
      timing_data.global_stop_count=0;
      /* CLK_TCK comes from the bits/time.h include file which is included
         in time.h which is included in djb_timing_b.h. It is supposed to
	 be the number of clock ticks per second (proc speed).
	 But it obviously isn't as its value is 1024. We'll try 1024^3
      timing_data.secs_per_clock = d_one/((double)CLK_TCK);
      */
      timing_data.secs_per_clock = d_one/((double)(1073741824L));
      timing_data.index = TIMING_MAX_ENTRIES;
      
      j = 0;
      for (i=0;i<TIMING_MAX_ENTRIES;i++) {
  	(timing_data.timings)[i].acc_time=(int64_t)0;
  	(timing_data.timings)[i].local_start_count=(int64_t)0;
  	(timing_data.timings)[i].local_stop_count=(int64_t)0;
	(timing_data.timings)[i].call_count=(int64_t)0;
  	(timing_data.timings)[i].tagcs = &(timing_data.timings_alpha[j]);
	(timing_data.timings)[i].indent = 0;
	(timing_data.timings_alpha)[j] = '\0';
  	j += 1 + TIMING_MAX_TAG_LEN;
      }
      TIMINGI_FP = fopen(timings_list_indx,"r");
      if (TIMINGI_FP == NULL) {
	fprintf(stderr,"Unable to open timings index file, %s\n",timings_list_indx);
	timing_data.alpha_present = 0;
      } else {
	timing_data.alpha_present = 1;
	for (i=DJB_TIMING_MIN;i<TIMING_MAX_ENTRIES;i++) {
	  fgets((char*)timingi_buf,TIMINGI_BUF_LEN,TIMINGI_FP);
	  if (feof(TIMINGI_FP)) break;
	  noents = sscanf((char*)timingi_buf,"%s %s %d,%d",
			  (char*)def_str,(char*)desc,&val,&indent);
	  if (noents == 4) {
	    if (strncmp((char*)def_str,"#define",7) == 0) {
	      val = val & DJB_TIMING_MAX;
	      /*
	      if (val >= TIMING_MAX_ENTRIES) {
		fprintf(stderr,
		    "Warning: from timing_init: Too many timing entries.\n");
		fprintf(stderr,
	   " Please increase TIMING_MAX_ENTRIES in djb_timing_b.h and recompile.\n");
		val = 0;
	      } else {
	      */
	      if (strlen((char*)desc) > TIMING_MAX_TAG_LEN) {
		strncpy((timing_data.timings[val]).tagcs,(char*)desc,
			TIMING_MAX_TAG_LEN);
		((timing_data.timings)[val].tagcs)[TIMING_MAX_TAG_LEN]='\0';
	      } else {
		strcpy((timing_data.timings)[val].tagcs,(char*)desc);
	      }
	      timing_data.timings[val].indent = indent;
              /*		
	      }
	      */
	    }
	  }
	}
	fclose(TIMINGI_FP);
      }
      /* 
  	Time the timing_start/timing_stop routines.
      */
      timing_data.start_stop_index = DJB_TIMING_START_STOP;
      timing_data.init_index = DJB_TIMING_INIT;
      timing_data.envelope_index = DJB_TIMING_ENVELOPE;
      scr_indent = 0;
      djb_timing_start(DJB_TIMING_START_STOP,scr_indent);
      for (i=0;i<TIMING_ITERATES;i++) {
  	djb_timing_start(DJB_TIMING_ENVELOPE,scr_indent);
  	djb_timing_stop(DJB_TIMING_ENVELOPE,scr_indent);
      }
      djb_timing_stop(DJB_TIMING_START_STOP,scr_indent);
      timing_data.ticks_per_envelopes = 
        timing_data.timings[DJB_TIMING_ENVELOPE].acc_time;
      /*
      fprintf(stdout,"ticks_per_envelope = %ld\n",
              timing_data.ticks_per_envelopes);
      fflush(stdout);
      */
      /* 
        do the mumble, to let the timer call to itc_clock finish.
      */
      stopt = itc_clock();
      mumble_cnt = 7;
      mumble1 = start;
      mumble2 = timing_data.ticks_per_envelopes;
      for (i=0;i<mumble_cnt;i++) {
	mumble2 = mumble1 + mumble2;
	mumble1 = (mumble2 << 23) | (mumble2 >> 41);
      }
      if (mumble1 == 42) {
	fprintf(stderr," mumble = 42\n");
	fflush(stderr);
      }
      timing_data.timings[DJB_TIMING_INIT].acc_time = stopt - start;
    } 
    /*
  }
    */
}






