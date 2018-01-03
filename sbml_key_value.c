/* sbml_key_value.c
  
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
*/
#include "boltzmann_structs.h"

#include "sbml_key_value.h"

int sbml_key_value(char *buffer, char *key, char *value, int max_key_len, 
		    int max_value_len) {

  /*
    scan a buffer for keyword "=" "string" triples returning
    full length allowing whitespace before keyword, between any of the
    three tokens and replacing all whitespace in string with underscores.
    Called by: parse_sbml
    Calls:     count_ws,count_nws,strlen,strncpy
  */
  char *line;
  int tl;
  int nb;

  int kl;
  int vl;

  int i;
  int ll;

  int eq_pos;
  int q_pos;

  line = buffer;
  tl = 0;
  ll   = strlen(buffer);
  nb = count_ws(line);
  tl = nb;
  line += nb; /* caution address arithmetic. */
  eq_pos = 0;
  /*
    Find the equals sign.
  */
  for (i=0;i<ll-tl;i++) {
    if (line[i] == '=') {
      eq_pos = i;
      break;
    }
  }
  /*
    Compute the key length allowing that there may or may not be white space
    between the key and the equal sign.
  */
  kl = count_nws(line);
  if (kl > eq_pos) {
    kl = eq_pos - 1;
  }
  if (kl > 0) {
    /*
      copy the key in to the key variable.
    */
    if (kl < max_key_len) {
      strncpy(key,line,kl);
      key[kl] = '\0';
      tl += eq_pos + 1;
      line += eq_pos + 1; /* caution address arithmetic. */
      /*
	Look for an opening double quote.
      */
      nb = count_ws(line);
      tl += nb;
      line += nb; /* caution address arithmetic. */
      q_pos = -1;
      for (i=0;i<(ll-tl);i++) {
	if (line[i] == '"') {
	  q_pos = i;
	  break;
	}
      }
      if (q_pos >= 0) {
	/*
	  found opening quote, skip vover it.
	*/
	line += q_pos + 1; /* caution address arithmetic. */
	tl += q_pos + 1;
	/*
	  Now find closing quote.
	*/
	q_pos = -1;
	for (i=0;i<ll-tl;i++) {
	  if (line[i] == '"') {
	    q_pos = i;
	    break;
	  }
	}
	if (q_pos > 0) {
	  /*
	    so value length in q_pos.
	  */
	  vl = q_pos;
	  if (vl < max_value_len) {
	    strncpy(value,(char*)&line[0],vl);
	    value[vl] = '\0';
	    tl += vl + 1;
	    /*
	      Now we want to strip trailing and leading blandks from
	      value and replace interior whitespace with underscores.
	    */
	    /*
	      remove trailing blanks.
	    */
	    for (i=vl-1;i>=0;i--) {
	      if (value[i] == ' ') {
		value[i] = '\0';
	      } else {
		break;
	      }
	    }
	    /*
	      Remove leading blanks.
	    */
	    vl = strlen(value);
	    nb = count_ws(value);
	    if (nb  > 0) {
	      for (i=0;i<(vl-nb);i++) {
		value[i] = value[i+nb];
	      }
	      vl = vl - nb;
	    }
	    /*
	      Replace interior blanks with underscores.
	    */
	    for (i=1;i<vl-1;i++) {
	      if (value[i] == ' ') {
		value[i] = '_';
	      }
	    }
	  } else {
	    /*
	      Value too long.
	    */
	    /*
	      fprintf(stderr,"parse_key_value: Error value too long, %s\n",
	              line);
	      fflush(stderr);
	    */
	    tl = -1;
	  }
	} else {
	  /*
	    Did not find close quote on value, or empty value.
	  */
	  /*
	    fprintf(stderr,"parse_key_value: Error missing closing quote on "
	    "key value pair, or value was empty string\n");
	    fflush(stderr);
	  */
	  tl = -1;
	}
      } else {
	/*
	  Did not find open quote on value.
	*/
	/*
	  fprintf(stderr,"parse_key_value: Error missing opening quote on "
	  "key value pair\n");
	  fflush(stderr);
	*/
	tl = -1;
      }
    } else {
      /*
	Key length is too long.
      */
      /*
	fprintf(stderr,"parse_key_value: Error keyword too long %s\n",
	        line);
	fflush(stderr);
      */
      tl = -1;
    }

  } else {
    /*
      Did not find an equal sign, error.
    */
    tl = -1;
    /*
      fprintf(stderr,"parse_key_value: Error missing equal sign, = \n");
      fflush(stderr);
    */
  }
  return (tl);
}

