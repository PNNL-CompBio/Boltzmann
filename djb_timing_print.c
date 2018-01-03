/* djb_timing_print.c
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

void djb_timing_print(FILE *fp) {
  /*
    Initialize the timings data structure.
  */
  double  d_zero, d_one, ooit, tperstartstop, tperenvelope, tinit;
  double  atime, rtime, ctime, etime;
  double  spc,ooitpspc;
  int64_t tperstartstops,tperenvelopes;
  int64_t start, stopt;
  int64_t startc,stopc;
  int64_t callc;
  char *tag;
  char label[64];
  int  ssi, initi, i, j;
  int  level, tab;
  char indent[2*TIMING_MAX_INDENT+1];
  /*
  fprintf (stdout, "djb_timing_print: &timing_data = %px\n",&timing_data);
  fflush(stdout);
  */
  start = itc_clock();
  d_one = 1.0;
  d_zero= 0.0;
  for (i=0;i<2*TIMING_MAX_INDENT+1;i++){
    indent[i] = ' ';
  }
  /*
    Need to check that timing_data is a valid pointer.
  */
  /*
  if (timing_data == NULL) {
    fprintf(stderr," Error: in timing_print: timing_data pointer is NULL.\n");
    fprintf(stderr,
	    " You need to include \"timing_init.h\" in your main program.\n");
  } else {
  */
    if (timing_data.initialized != TIMING_INITIALIZED) {
      fprintf(stderr,
	      " Error: in timing_print: timing_data was not initialized.\n");
      fprintf(stderr,
	    " You need to include \"timing_init.h\" in your main program.\n");
    } else {
      ooit = d_one/((double)TIMING_ITERATES);
      spc  = timing_data.secs_per_clock;
      ooitpspc = ooit * spc;
      ssi = timing_data.start_stop_index;
      initi = timing_data.init_index;
      tperstartstops = (timing_data.timings)[DJB_TIMING_START_STOP].acc_time;
      tperenvelopes  = tperstartstops - 
	(timing_data.timings)[DJB_TIMING_ENVELOPE].acc_time;
      tperstartstop = (double)(tperstartstops) * ooitpspc;
      tperenvelope  = (double)(tperenvelopes) * ooitpspc;
      tinit         = (double)(timing_data.timings)[DJB_TIMING_INIT].acc_time*spc;

      /*
        Print timer start, stop and init data.
      */
      if (fp) {
	fprintf(fp," Seconds per clock    = %12.5e\n",spc);
	fprintf(fp," Seconds per envelope = %12.5e\n",tperenvelope);
	fprintf(fp,
		"Timing Tag                            Raw_Time        CallC "
		"  Start/StopC  Act_Time\n");
	rtime = (double)(timing_data.timings)[DJB_TIMING_START_STOP].acc_time*spc;
	fprintf(fp,"%-35s %12.5e %10ld %10ld %12.5e\n",
		"   DJB_TIME_START_STOP",rtime,(long)1,
		(long)TIMING_ITERATES,tperstartstop);
	fprintf(fp,"%-35s %12.5e %10ld %10ld %12.5e\n",
		"   DJB_TIME_INIT",tinit,(long)1,
	      (long)TIMING_ITERATES,tinit);
	for (i=DJB_TIMING_MIN; i<=DJB_TIMING_MAX;i++) {
	  tag    = (timing_data.timings)[i].tagcs;
	  callc = (timing_data.timings)[i].call_count;
	  startc = (timing_data.timings)[i].local_start_count;
	  stopc = (timing_data.timings)[i].local_stop_count;
	  rtime = (double)((timing_data.timings)[i].acc_time)*spc;
	  etime = (((double)callc)*((double)tperenvelopes))*ooitpspc;
	  ctime = (((double)startc)*((double)tperstartstops))*ooitpspc;
	  atime = rtime - ctime - etime;
	  level = (timing_data.timings)[i].indent;
	  tab = level+level;
	  if (tab > 20) {
	    tab = 20;
	  }
	  indent[tab] = '\0';
	  strcpy(label,indent);
	  indent[tab] = ' ';
	  sprintf((char *)&label[tab],"%2d ",level);
	  strcpy((char*)&label[tab+3],tag);
	  /*	if ((strlen(tag) != 0) || (atime > 0.0)) {}*/
	  if ((strlen(tag) != 0) && (atime > 0.0)) {
	    if (startc != stopc) {
	      fprintf(fp,
		"unbalanced start/stop inside %s startc = %lld stopc = %lld\n",
		      tag,startc,stopc);
	    }
	    fprintf(fp,"%-35s %12.5e %10lld %10lld %12.5e\n",
		    label,rtime,callc,startc,atime);

	  }
	}
	stopt = itc_clock();
	rtime = (double)(stopt-start)*spc;
	fprintf(fp,"%-35s %12.5e %10lld %10lld %12.5e\n",
		"   DJB_TIMING_PRINT",rtime,(int64_t)0,(int64_t)0,rtime);
      } /* end if (fp) */
    } /* end else Timing was initialized. */
  /*
  } // end else Non null timing_data pointer 
  */

}
