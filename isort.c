#include "boltzmann_structs.h"
#include "imerge.h"
#include "isort.h"
void isort(int n,int *list, int *scratch) {
  /*
    Sort a list of n int's using scratch space of length n, and
    overwriting the input list with the sorted result.
    Called by: iluvf
    Calls:     imerge, memcpy
  */
  int64_t move_size;
  int64_t l_4;
  int *list1;
  int *list2;
  int *itemp;

  int l1;
  int l2;

  int step;
  int j;

  int next_step;
  int ln;

  int k;
  int padi;

  list1 = list;
  list2 = scratch;
  l_4   = (int64_t)4;
  
  if (n > 1) {
    for (step = 1; step < n; step += step) {
      next_step = step + step;
      for (j=0;j<(n-step); j+= next_step) {
	l1 = step;
	l2 = n - j - step;
	if (l2 > step) {
	  l2 = step;
	}
	imerge((int*)&list1[j],(int*)&list1[j+step],(int*)&list2[j],l1,l2);
      } /* end for j */
      ln = n & (next_step - 1);
      if (ln <= step) {
	move_size = ln * l_4;
	k = n - ln;
	if (move_size > 0) {
	  memcpy((void*)&list2[k],(void*)&list1[k],move_size);
	}
      }
      itemp = list1;
      list1 = list2;
      list2 = itemp;
    } /* end for step */
  } /* end if (n > 1) */
  if (list1 != list) {
    move_size = n * l_4;
    memcpy((void*)list,(void*)list1,move_size);
  }
}
