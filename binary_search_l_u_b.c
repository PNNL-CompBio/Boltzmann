/* binary_search_l_u_b.c
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

#include "binary_search_l_u_b.h"
int binary_search_l_u_b (double *d, double v, int n) {
  /*
    find smallest j such that d[j] >= v;
    assumes d[n-1] > v, for (0<=j<n);

    Called by : boltzmann_run
    Calls"
  */
  int low;
  int mid;
  int high;
  int result;

  if (v <= d[0]) {
    result = 0;
  } else {
    if (v > d[n-2]) {
      result = n-1;
    } else {
      low = 0;
      high = n-1;
      mid = (low+high)>>1;
      while (mid != low) {
	if (d[mid] >= v) {
	  high = mid;
	} else {
	  low = mid;
	}
	mid = (low+high)>>1;
      }
      result = high;
    }
  }
  return (result);
}
