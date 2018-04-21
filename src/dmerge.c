#include "boltzmann_structs.h"
#include "dmerge.h"
void dmerge(double *list1, double *list2, double *list3, int l1, int l2) {
  /*
    Merge two sorted arrays of int's, list1 and list2, into a sorted single
    list, list 3.
    Called by: isort
    Calls:     memcpy
  */
  int64_t move_size;
  int64_t l_8;
  double *dp1;
  double *dp2;
  double *dp3;
  int j1;
  int j2;
  int j3;
  int n;
  l_8 = (int64_t)8;
  n = l1 + l2;
  j1 = 0;
  j2 = 0;
  j3 = 0;
  dp1 = list1;
  dp2 = list2;
  dp3 = list3;
  for (j3 = 0; j3 < n; j3++ ) {
    if (*dp1 <= *dp2) {
      *dp3 = *dp1;
      dp1 += 1; /* Caution address arithmetic here. */
      dp3 += 1; /* Caution address arithmetic here. */
      j1++;
      if (j1 == l1) {
	move_size = (l2 - j2) * l_8;
	if (move_size > 0) {
	  memcpy(dp3,dp2,move_size);
	}
	break;
      }
    } else {
      *dp3 = *dp2;
      dp2 += 1; /* Caution address arithmetic here. */
      dp3 += 1; /* Caution address arithmetic here. */
      j2 ++;
      if (j2 == l2) {
	move_size = (l1 - j1) * l_8;
	if (move_size > 0) {
	  memcpy(dp3,dp1,move_size);
	}
	break;
      }
    }
  } /* end for (j3...) */
}
