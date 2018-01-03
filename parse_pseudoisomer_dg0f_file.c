/* parse_pseudoisomer_dg0f_file.c
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
/*
 * parse_pseudoisomer_dg0f_file.c
 *
 *  Created on: Apr 17, 2013
 *      Author:  Dennis G. Thomas
 *  Adapted by Doug Baxter May 24, 2013
*/
#include "boltzmann_structs.h"

#include "upcase.h"
#include "blank_to_dash.h"
#include "sharp_pos.h"

#include "parse_pseudoisomer_dg0f_file.h"
int parse_pseudoisomer_dg0f_file(struct pseudoisomer_struct *pseudoisomers,
				 char *pseudoisomer_strings,
				 char *pseudoisomer_file,
				 int64_t  num_cpd,
				 int64_t  align_len) {
  /*
    Parses the energy of formation file building the pseudoisomers
    structure (an array of pseudoisomer_struct's) and the 
    pseudoiosomer_strings array.
    Called by: compute_standard_energies
    Calls:     sharp_pos, fopen, fgets, strcpy, fprintf, fflush, sscanf,
               fclose
  */
  struct pseudoisomer_struct *pseudoisomer;
  int64_t s_pos;
  int64_t padded_len;
  int64_t align_mask;
  int64_t cpd_ct;
  double z;
  double dg0f;
  /*
  char *pseudoisomer_strings;
  char *pseudoisomer_file;
  */
  char *buffer;
  char *bptr;
  char* next_string;
  char buffer_v[1024];

  int success;
  int buff_len;

  int i;
  int tok_len;
  
  int nh;
  int priority_level;

  int ns;
  int padi;

  FILE *psi_fp;

  success = 1;
  /*
  pseudoisomers        = state->pseudoisomers;
  pseudoisomer_strings = state->pseudoisomer_strings;
  pseudoisomer_file    = state->pseudoisomer_dg0f_file;
  num_cpd              = state->num_pseudoisomers;
  align_len            = state->align_len;
  */
  align_mask           = align_len - (int64_t)1;
  s_pos                = (int64_t)0;
  next_string = (char *)&pseudoisomer_strings[s_pos];
  buff_len = 1024;
  buffer   = (char *)&buffer_v[0];
  
  /* 
    open input file for reading 
  */
  psi_fp = fopen(pseudoisomer_file,"r");

  if(psi_fp  == NULL){
    success = 0;
    fprintf(stderr,"parse_pseudoisomer_dg0f_file: pseudoisomer_dg0f file, %s,"
	    " not open.\n",
	    pseudoisomer_file);
    fflush(stderr);
  }
  /*
    chew up the header line.
  */
  cpd_ct = (int64_t)0;
  if (success) {
    bptr = fgets(buffer,buff_len,psi_fp);
    if (bptr) {
      cpd_ct = (int64_t)0;
      pseudoisomer = pseudoisomers;
      for (i=0;i<num_cpd;i++) {
	bptr = fgets(buffer,buff_len,psi_fp);
	if (bptr == NULL) {
	  fprintf(stderr,"parse_pseudoisomer_dg0f_file: Error premature "
		  "end of file \n");
	  fflush(stderr);
	  success = 0;
	  break;
	}
	/* 
	   get the json cpd name 
	*/	  
	tok_len = sharp_pos(bptr);
	if (tok_len < 1) {
	  fprintf(stderr,"parse_pseudoisomer_dg0f_file: Error reading "
		  "json cpd name in line %d, line was\n%s\n",i,buffer);
	  fflush(stderr);
	  continue;
	}   
	bptr[tok_len] = '\0';
	tok_len += 1;
	strcpy(next_string,bptr);
	upcase(tok_len,next_string,next_string);
	blank_to_dash(tok_len,next_string,next_string);
	padded_len = ((align_len - (tok_len & align_mask)) & align_mask) + 
  	             tok_len;
	pseudoisomer->json_cpd_name = s_pos;
	next_string += padded_len; /* Caution address arithmetic. */
	bptr  += tok_len;          /* Caution address arithmetic. */
	s_pos += padded_len;
	
	/*
	  Get the kegg ID
	*/
	tok_len = sharp_pos(bptr);
	if (tok_len < 1) {
	  fprintf(stderr,"parse_pseudoisomer_dg0f_file: Error reading "
		  "kegg id in line %d, line was\n%s\n",i,buffer);
	  fflush(stderr);
	  continue;
	}   
	bptr[tok_len] = '\0';
	tok_len += 1;
	strcpy(next_string,bptr);
	padded_len = ((align_len - (tok_len & align_mask)) & align_mask) + 
  	             tok_len;
	pseudoisomer->kegg_id = s_pos;
	next_string += padded_len; /* Caution address arithmetic. */
	bptr  += tok_len;          /* Caution address arithmetic. */
	s_pos += padded_len;
	/* 
	  get the priority level 
	*/
	tok_len = sharp_pos(bptr);
	if (tok_len < 1) {
	  fprintf(stderr,"parse_pseudoisomer_dg0f_file: Error reading "
		  "priority_level in line %d, line was\n%s\n",i,buffer);
	  fflush(stderr);
	  continue;
	}   
	bptr[tok_len] = '\0';
	tok_len += 1;
	ns = sscanf(bptr,"%d",&priority_level);
	if (ns == 1) {
	  pseudoisomer->priority_level = priority_level;
	} else {
	  fprintf(stderr,"parse_pseudoisomer_dg0f_file: Error reading "
		  "priority level in line %d, line was\n%s\n",i,buffer);
	  fflush(stderr);
	  continue;
	}   
	bptr  += tok_len;       /* Caution address arithmetic. */
	/* 
	  get the source name 
	*/
	tok_len = sharp_pos(bptr);
	if (tok_len < 1) {
	  fprintf(stderr,"parse_pseudoisomer_dg0f_file: Error reading "
		  "source name in line %d, line was\n%s\n",i,buffer);
	  fflush(stderr);
	  continue;
	}   
	bptr[tok_len] = '\0';
	tok_len += 1;
	strcpy(next_string,bptr);
	padded_len = ((align_len - (tok_len & align_mask)) & align_mask) + 
  	             tok_len;
	pseudoisomer->source = s_pos;
	next_string += padded_len; /* Caution address arithmetic. */
	bptr  += tok_len;          /* Caution address arithmetic. */
	s_pos += padded_len;
	/* 
	  get the value of nH 
	*/
	tok_len = sharp_pos(bptr);
	if (tok_len < 1) {
	  fprintf(stderr,"parse_pseudoisomer_dg0f_file: Error reading nh "
		  "in line %d, line was\n%s\n",i,buffer);
	  fflush(stderr);
	  continue;
	}   
	bptr[tok_len] = '\0';
	tok_len += 1;
	ns = sscanf(bptr,"%d",&nh);
	if (ns == 1) {
	  pseudoisomer->nh = nh;
	} else {
	  fprintf(stderr,"parse_pseudoisomer_dg0f_file: Error reading nh "
		  "in line %d, line was\n%s\n",i,buffer);
	  fflush(stderr);
	  fflush(stderr);
	  continue;
	}   
	bptr  += tok_len;       /* Caution address arithmetic. */

	/* 
	  get the value of z the ionic strength.
	*/
	tok_len = sharp_pos(bptr);
	if (tok_len < 1) {
	  fprintf(stderr,"parse_pseudoisomer_dg0f_file: Error reading "
		  "ionic strength, z, in line %d, line was\n%s\n",i,buffer);
	  fflush(stderr);
	  continue;
	}   
	bptr[tok_len] = '\0';
	tok_len += 1;
	ns = sscanf(bptr,"%le",&z);
	if (ns == 1) {
	  pseudoisomer->z = z;
	} else {
	  fprintf(stderr,"parse_pseudoisomer_dg0f_file: Error reading "
		  "ionic strength, z, in line %d, line was\n%s\n",i,buffer);
	  fflush(stderr);
	  continue;
	}   
	bptr  += tok_len;       /* Caution address arithmetic. */
	/* 
	  get the value of delta G0_f (kJ/mol) 
	*/
	tok_len = sharp_pos(bptr);
	if (tok_len < 1) {
	  fprintf(stderr,"parse_pseudoisomer_dg0f_file: Error reading "
		  "dg0f in line %d, line was\n%s\n",i,buffer);
	  fflush(stderr);
	  continue;
	}   
	bptr[tok_len] = '\0';
	tok_len += 1;
	ns = sscanf(bptr,"%le",&dg0f);
	if (ns == 1) {
	  pseudoisomer->dg0_f = dg0f;
	} else {
	  fprintf(stderr,"parse_pseudoisomer_dg0f_file: Error reading "
		  "dg0f in line %d, line was\n%s\n",i,buffer);
	  fflush(stderr);
	  continue;
	}   
	bptr  += tok_len;       /* Caution address arithmetic. */
	/*
	  Get the reference string.
	*/
	tok_len = strlen(bptr);
	if (tok_len < 1) {
	  fprintf(stderr,"parse_pseudoisomer_dg0f_file: Error reading "
		  "reference string in line %d, line was\n%s\n",i,buffer);
	  fflush(stderr);
	  continue;
	}   
	bptr[tok_len-1] = '\0';
	strcpy(next_string,bptr);
	padded_len = ((align_len - (tok_len & align_mask)) & align_mask) + 
  	             tok_len;
	pseudoisomer->ref = s_pos;
	next_string += padded_len; /* Caution address arithmetic. */
	s_pos += padded_len;

	pseudoisomer += 1;         /* Caution address arithmetic. */
	cpd_ct += 1;
      } /* end for(i...) */
    } /* end if (bptr) */
    fclose(psi_fp);
  } /* end if (success) */
  if (cpd_ct != num_cpd) {
    success = 0;
  }
  return (success);
}


