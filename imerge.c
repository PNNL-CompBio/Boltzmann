#include "boltzmann_structs.h"
#include "imerge.h"
void imerge(int *list1, int *list2, int *list3, int l1, int l2) {
  /*
    Merge two sorted arrays of int's, list1 and list2, into a sorted single
    list, list 3.
    Called by: isort
    Calls:     memcpy
  */
  int64_t move_size;
  int64_t l_4;
  int *ip1;
  int *ip2;
  int *ip3;
  int j1;
  int j2;
  int j3;
  int n;
  l_4 = (int64_t)4;
  n = l1 + l2;
  j1 = 0;
  j2 = 0;
  j3 = 0;
  ip1 = list1;
  ip2 = list2;
  ip3 = list3;
  for (j3 = 0; j3 < n; j3++ ) {
    if (*ip1 <= *ip2) {
      *ip3 = *ip1;
      ip1 += 1; /* Caution address arithmetic here. */
      ip3 += 1; /* Caution address arithmetic here. */
      j1++;
      if (j1 == l1) {
	move_size = (l2 - j2) * l_4;
	if (move_size > 0) {
	  memcpy(ip3,ip2,move_size);
	}
	break;
      }
    } else {
      *ip3 = *ip2;
      ip2 += 1; /* Caution address arithmetic here. */
      ip3 += 1; /* Caution address arithmetic here. */
      j2 ++;
      if (j2 == l2) {
	move_size = (l1 - j1) * l_4;
	if (move_size > 0) {
	  memcpy(ip3,ip1,move_size);
	}
	break;
      }
    }
  } /* end for (j3...) */
}
