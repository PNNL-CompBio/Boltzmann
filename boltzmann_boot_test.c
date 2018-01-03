/* boltzmann_boot_test.c
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
/*
#define DBG_BOLTZMANN_BOOT  
*/
#include "boltzmann_boot.h"
#include "boltzmann_boot_check.h"

int main(int argc, char **argv)
{
  struct super_state_struct *super_statep;
  char    *param_file_name;
  int     success;
  int     retval;
  FILE    *lfp;
  if (argc > 1) {
    param_file_name = argv[1];
    success =  boltzmann_boot(param_file_name,&super_statep);

  } else {
    fprintf(stderr,"boltzmann_boot_test: Error no param file specified as argument\n");
    fflush(stderr);
    param_file_name = NULL;
    success = 0;
  }
  if (success) {
    /*
      Print out each of the reaction files and their initial concentrations.
    */
    lfp = fopen("boot_check.out","w");
    if (lfp) {
      success = boltzmann_boot_check(super_statep,lfp);
      fclose(lfp);
    } else {
      fprintf(stderr,"boltzmann_boot_test: Error could not open boot_check.out\n");
      fflush(stderr);
    }
  }
  if (success) {
    retval = 0;
  } else {
    retval = 1;
  }
  return(retval);
}
