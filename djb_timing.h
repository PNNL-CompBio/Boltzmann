/* djb_timing.h
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
  This file is included in the timingi.h file
  as defines te timing macros for using routines in libdjb_timing.a
*/
#ifndef _djb_timing_h
#define _djb_timing_h

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
extern int64_t itc_clock(void);

#if defined(TIMING_ON) && !defined(_lint)

extern void djb_timing_init(char *timings_list_indx);

extern void djb_timing_start(int index, int indent);


extern void djb_timing_stop(int index, int indent);


extern void djb_timing_print(FILE *fp);

extern void djb_timing_reset(void);

#define TIMING_INIT(a) djb_timing_init(a);
#define TIMING_START(a) djb_timing_start(a);
#define TIMING_STOP(a)  djb_timing_stop(a);
#define TIMING_PRINT(a) djb_timing_print(a);
#define TIMING_RESET    djb_timing_reset();

#else
#define TIMING_INIT(a) {};
#define TIMING_START(a) {};
#define TIMING_STOP(a)  {};
#define TIMING_PRINT(a) {};
#define TIMING_RESET {};
#endif
#include "timingi.h"
#ifdef	__cplusplus
}
#endif
#endif
