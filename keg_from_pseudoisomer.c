/* keg_from_pseudoisomer.c
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
#include "check_for_ws.h"
#include "boltzmannize_json_id.h"

int main(int argc, char **argv) {
  /*
    From the pseudoisomer_dg0f.txt file produce of file of unique json_id's
    sorted by keg_id.
    produces file unique_json_ids_sorted_by_keg.txt
    Calls: fopen, fprintf, fflush, fgets, fclose,
           count_ws,
	   count_nlb,
	   check_for_ws,
	   boltzmannize_json_id
  */
  char buff_c[67108864];
  char *buffer;
  char *line;
  char *json_id;
  char *keg_id;
  int  buff_len;
  int  ns;
  int  nc;
  int  json_id_has_ws;
  FILE *infp;
  FILE *json_w_ws_fp;
  FILE *unsorted_ids_fp;
  infp = fopen("pseudoisomer_dg0f.txt","r");
  json_w_ws_fp = fopen("unsorted_json_ids_with_spaces.txt","w");
  unsorted_ids_fp = fopen("unsorted_ids.txt","w");
  if (infp == NULL) {
    fprintf(stderr,"keg from pseudoisomer could not open pseudoisomer_dg0f.txt\n");
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
      json_id = line;
      nc = count_nlb(line);
      line[nc] = '\0';
      line += (nc + 1); /* Caution address arithmetic */
      /*
	Skip white_space.
      */
      ns = count_ws(line);
      line += ns; /* Caution address arithmetic */
      keg_id = line;
      nc = count_nlb(line);
      line[nc] ='\0';

      json_id_has_ws = check_for_ws(line);
      if (json_id_has_ws) {
	fprintf(json_w_ws_fp,"%s\n",json_id);
      }
      /*
	Uppercase and replace white space with underscores.
      */
      boltzmannize_json_id(json_id);
      fprintf(unsorted_ids_fp,"%s\t%s\n",keg_id,json_id);
      line = fgets(buffer,buff_len,infp);
    } /* end while() */
    fclose(json_w_ws_fp);
    fclose(unsorted_ids_fp);
    system("sort -k 1 -u unsorted_ids.txt > unique_json_ids_sorted_by_keg.txt");
    system("sort -k 1 -u unsorted_json_ids_with_spaces.txt > json_ids_with_spaces.txt");
  }
  return(0);
}
