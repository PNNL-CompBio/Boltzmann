#include "boltzmann_structs.h"
#include "dcrsng_mag_merge.h"
void dcrsng_mag_merge(double *list1, int *index1,
		      double *list2, int *index2,
		      double *mlist, int *mindex,
		      int l1,
		      int l2) {
  /*
    Merge two lists of doubles in decreasing order of magnitude into a
    single list in decreasing order of magnitued with index fields 
    tagging along.
    Called by: dcrsng_mag_sort
    Calls:     memcpy
  */
  double *p1;
  double *p2;
  double *p3;
  int64_t l_8;
  int *ip1;
  int *ip2;
  int *ip3;
  int64_t move_size;
  
  int j1;
  int j2;
  int j3;
  int n;

  l_8 = (int64_t)8;
  j1 = 0;
  j2 = 0;
  j3 = 0;
  n = l1 + l2;
  p1 = list1;
  p2 = list2;
  p3 = mlist;
  ip1 = index1;
  ip2 = index2;
  ip3 = mindex;
  
  for (j3=0; j3 < n; j3++) {
    if (fabs(*p1) >= fabs(*p2)) {
      *p3  = *p1;
      *ip3 = *ip1;
      j1++;
      p1  += 1; /* Caution address arithmetic here. */
      ip1 += 1; /* Caution address arithmetic here. */
      p3  += 1; /* Caution address arithmetic here. */
      ip3 += 1; /* Caution address arithmetic here. */
      if (j1 == l1) {
	move_size = (l2-j2) * l_8;
	if (move_size > 0) {
	  memcpy(p3,p2,move_size);
	  move_size = move_size >> 1;
	  memcpy(ip3,ip2,move_size);
	}
	break;
      } 
    } else {
      *p3  = *p2;
      *ip3 = *ip2;
      j2++;
      p2  += 1; /* Caution address arithmetic here. */
      ip2 += 1; /* Caution address arithmetic here. */
      p3  += 1; /* Caution address arithmetic here. */
      ip3 += 1; /* Caution address arithmetic here. */
      if (j2 == l2) {
	move_size = (l1-j1) * l_8;
	if (move_size > 0) {
	  memcpy(p3,p1,move_size);
	  move_size = move_size >> 1;
	  memcpy(ip3,ip1,move_size);
	}
	break;
      } 
    } /* end else p1 had smaller magnitued */
  } /* end for(j3 ...) */
}
