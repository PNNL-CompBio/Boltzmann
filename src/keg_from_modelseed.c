/* keg_from_modelseed.c
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
#include "count_nlb.h"
#include "count_ntb.h"
#include "count_nor.h"
#include "check_for_ws.h"
#include "boltzmannize_json_id.h"

int main(int argc, char **argv) {
  /*
    From the modelseeed_compounds_db.txt file produce of file of unique 
    modelseed_id's sorted by keg_id.
    produces file unique_modelseed_ids_sorted_by_keg.txt
    and modelseed_multiple_keg_id.txt
    Calls: fopen, fprintf, fflush, fgets, fclose,
           count_ws,
	   count_ntb,
	   count_nor,
	   check_for_ws,
	   boltzmannize_json_id
  */
  char buff_c[67108864];
  char *buffer;
  char *line;
  char *modelseed_id;
  char *abbrev;
  char *names;
  char *name;
  char *keg_id;
  int  buff_len;
  int  ns;
  int  nc;
  int  nn;
  int  modelseed_id_has_ws;
  FILE *infp;
  FILE *modelseed_w_ws_fp;
  FILE *unsorted_ids_fp;
  FILE *multiple_keg_fp;
  infp = fopen("modelseed_compounds_db.txt","r");
  modelseed_w_ws_fp = fopen("unsorted_modelseed_ids_with_spaces.txt","w");
  unsorted_ids_fp = fopen("unsorted_modelseed_ids.txt","w");
  multiple_keg_fp = fopen("modelseed_multiple_keg_id.txt","w");
  if (infp == NULL) {
    fprintf(stderr,"keg from modelseed could not open modelseed_compounds_db.txt\n");
    fflush(stderr);
  } else {
    buff_len = 67108864;
    buffer = (char *)&buff_c[0];
    /*
      Chew up the title line.
    */
    line = fgets(buffer,buff_len,infp);
    line = fgets(buffer,buff_len,infp);
    while (!feof(infp)) {
      /*
	Skip initial white_space.
      */
      ns = count_ws(line);
      line += ns; /* Caution address arithmetic */
      /*
	Skip the database name.
      */
      nc = count_ntb(line);
      line += (nc + 1);
      modelseed_id = line;
      nc = count_ntb(line);
      line[nc] = '\0';
      line += (nc + 1); /* Caution address arithmetic */
      /*
	Get the abbreviation
      */
      abbrev = line;
      nc = count_ntb(line);
      line[nc] = '\0';
      line += (nc + 1); /* Caution address arithmetic */

      names = line;
      nc = count_ntb(line);
      line[nc] = '\0';
      line += (nc + 1); /* Caution address arithmetic */
      
      keg_id = line;
      nc = count_ntb(line);
      line[nc] ='\0';
      /*
	Check for multiple keg ids.
      */
      nn = count_nor(line);
      if (nn != nc) {
	fprintf(multiple_keg_fp,"%s\t%s\n",modelseed_id,keg_id);
	line[nn] = '\0';
      }
      line += (nc+1);
      

      modelseed_id_has_ws = check_for_ws(modelseed_id);
      if (modelseed_id_has_ws) {
	fprintf(modelseed_w_ws_fp,"%s\n",modelseed_id);
      }
      /*
	Uppercase and replace white space with underscores.
      */
      boltzmannize_json_id(modelseed_id);
      fprintf(unsorted_ids_fp,"%s\t%s\n",keg_id,modelseed_id);
      modelseed_id_has_ws = check_for_ws(abbrev);
      if (modelseed_id_has_ws) {
	fprintf(modelseed_w_ws_fp,"%s\n",abbrev);
      }
      boltzmannize_json_id(abbrev);
      fprintf(unsorted_ids_fp,"%s\t%s\n",keg_id,abbrev);
      /*
	Now need to loop over names.
      */
      name  = names;
      int cl;
      cl    = strlen(name);
      while (cl > 0) {
	nc = count_nor(name);
	name[nc] = '\0';
	cl -= (nc+1);
	modelseed_id_has_ws = check_for_ws(name);
	if (modelseed_id_has_ws) {
	  fprintf(modelseed_w_ws_fp,"%s\n",name);
	}
	boltzmannize_json_id(name);
	fprintf(unsorted_ids_fp,"%s\t%s\n",keg_id,name);
	name += (nc+1);
      }
      line = fgets(buffer,buff_len,infp);
    } /* end while() */
    fclose(modelseed_w_ws_fp);
    fclose(unsorted_ids_fp);
    fclose(multiple_keg_fp);
    system("sort -k 1 -u unsorted_modelseed_ids.txt > unique_modelseed_ids_sorted_by_keg.txt");
    system("sort -k 1 -u unsorted_modelseed_ids_with_spaces.txt > modelseed_ids_with_spaces.txt");
  }
  return(0);
}
