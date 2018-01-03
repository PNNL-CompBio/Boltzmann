/* modelseed_2_json.c
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
#include "boltzmann_structs.h"

#include "count_ws.h"
#include "count_nws.h"
#include "count_ntb.h"

int main(int argc, char **argv) {
  /*
    This routine creates the "modelseed_2_json.srt" file and the 
    "modelseed_json_mismatches" file from the 
    "kegg_2_modelseed.srt" and 
    "kegg_2_json.srt" files.
    It produces scratch file "modelseed_2_json.dat" (the unsorted file).
    
  */
  FILE *msfp;
  FILE *jsfp;
  FILE *m2jfp;
  FILE *mmfp;
  char mline_c[1048576];
  char jline_c[1048576];
  char *mline;
  char *jline;
  char *jp;
  char *mp;
  char *j_kegg_id;
  char *m_kegg_id;
  char *json_id;
  char *modelseed_id;
  int  line_len;
  int  nc;
  int  ns;
  line_len = 1048576;
  mline = (char *)&mline_c[0];
  jline = (char *)&jline_c[0];
  msfp = fopen("kegg_2_modelseed.srt","r");
  if (msfp != NULL) {
    jsfp = fopen("kegg_2_json.srt","r");
    if (jsfp != NULL) {
      m2jfp = fopen("modelseed_2_json.dat","w");
      if (m2jfp != NULL) {
	mmfp = fopen("modelseed_json_mismatches","w");
	if (mmfp != NULL) {
	  jp = fgets(jline,line_len,jsfp);
	  /*
	    Skip leading whitespace.
	  */
	  ns = count_ws(jp);
	  jp += ns; /* Caution address arithmetic */
	  j_kegg_id = jp;
	  nc = count_ntb(jp);
	  jp[nc] = '\0';
	  jp += (nc+1); /* Caution address arithmetic */
	  json_id = jp;
	  nc = count_nws(jp);
	  jp[nc] = '\0';
	  mp = fgets(mline,line_len,msfp);
	  while (!feof(msfp)) {
	    /*
	      Skip leading whitespace.
	    */
	    ns = count_ws(mp);
	    mp += ns; /* Caution address arithmetic */
	    m_kegg_id = mp;
	    nc = count_ntb(mp);
	    mp[nc] = '\0';
	    mp += (nc+1); /* Caution address arithmetic */
	    modelseed_id = mp;
	    nc = count_nws(mp);
	    mp[nc] = '\0';
	    while ((strcmp(j_kegg_id,m_kegg_id) < 0) && !feof(jsfp)) {
	      jp = fgets(jline,line_len,jsfp);
	      if (!feof(jsfp)) {
		ns = count_ws(jp);
		jp += ns; /* Caution address arithmetic */
		j_kegg_id = jp;
		nc = count_ntb(jp);
		jp[nc] = '\0';
		jp += (nc+1); /* Caution address arithmetic */
		json_id = jp;
		nc = count_nws(jp);
		jp[nc] = '\0';
	      }
	    }
	    if (strcmp(j_kegg_id,m_kegg_id) == 0) {
	      fprintf(m2jfp,"%s\t%s\n",modelseed_id,json_id);
	    } else {
	      fprintf(mmfp,"%s\t%s\n",m_kegg_id,modelseed_id);
	    }
	    mp = fgets(mline,line_len,msfp);
	  }
	} else {
	  fprintf(stderr,"unable to open modelseed_json_mismatches\n");
	}
      } else {
	fprintf(stderr,"unable to open modelseed_2_json.dat\n");
      }
    } else {
      fprintf(stderr,"unable to open kegg_2_json.srt\n");
    }
  } else {
    fprintf(stderr,"unable to open kegg_2_modelseed.srt\n");
  }
  fclose(msfp);
  fclose(jsfp);
  fclose(m2jfp);
  fclose(mmfp);
  system("sort -k 1 -u modelseed_2_json.dat > modelseed_2_json.srt");  
  return(0);
}
