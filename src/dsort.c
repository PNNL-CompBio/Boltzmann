#include "boltzmann_structs.h"
#include "dsort_pairs.h"
#include "dsort_pairs_in_place.h"
#include "dmerge.h"
#include "dsort.h"
void dsort(int n,double *list, double *scratch) {
  /*
    Sort a list of n doubles using scratch space of length n, and
    overwriting the input list with the sorted result.
    Called by: lr11_approximate_delta_concs, lr12_approximate_delta_concs
    Calls:     dmerge, memcpy, dsort_pairs, dsort_pairs_in_place
  */
  int64_t move_size;
  int64_t l_8;
  double *list1;
  double *list2;
  double *dtemp;

  int l1;
  int l2;

  int step;
  int j;

  int next_step;
  int ln;

  int k;
  int in_place_pairwise;

  list1 = list;
  list2 = scratch;
  l_8   = (int64_t)8;
  
  /*
    We want to avoid the last copy if possible.
  */
  in_place_pairwise = 0;
  if (n > 1) {
    for (step = 1;step < n; step += step) {
      in_place_pairwise = 1 - in_place_pairwise;
    }
    if (in_place_pairwise) {
       dsort_pairs_in_place(n,list);
       list1 = list;
       list2 = scratch;
    } else {
       dsort_pairs(n,list,scratch);
       list1 = scratch;
       list2 = list;
    }
  }
  /*
    NB the following loop only executes if n > 2.
  */
  for (step = 2; step < n; step += step) {
    next_step = step + step;
    for (j=0;j<(n-step); j+= next_step) {
      l1 = step;
      l2 = n - j - step;
      if (l2 > step) {
	l2 = step;
      }
      dmerge((double*)&list1[j],(double*)&list1[j+step],(double*)&list2[j],l1,l2);
    } /* end for j */
    ln = n & (next_step - 1);
    if (ln <= step) {
      move_size = ln * l_8;
      k = n - ln;
      if (move_size > 0) {
	memcpy((void*)&list2[k],(void*)&list1[k],move_size);
      }
    }
    dtemp = list1;
    list1 = list2;
    list2 = dtemp;
  } /* end for step */
  /*
  if (list1 != list) {
    move_size = n * l_8;
    memcpy((void*)list,(void*)list1,move_size);
  }
  */
}
