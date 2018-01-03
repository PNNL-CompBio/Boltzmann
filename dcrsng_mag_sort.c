#include "boltzmann_structs.h"
#include "dcrsng_mag_merge.h"
#include "dcrsng_mag_sort.h"
void dcrsng_mag_sort(int n,
		     double *row,
		     int *index,
		     double *srow,
		     int *sindex) {
  /*
    Sort an array of doubles, row into decreasing order by magnitude, with
    the index array elements tagging along.
    srow and sindex are scratch arrays of length n.
    Called by: iluvf
    Calls:     dcrsng_mag_merge,memcpy
  */
  double *trow;
  double *tsrow;
  double *dtemp;
  int    *tindex;
  int    *tsindex;
  int    *itemp;
  int64_t move_size;
  int64_t l_8;
  int64_t l_4;

  int l1;
  int l2;
  
  int ln;
  int j;

  int k;
  int padi;

  int step;
  int next_step;

  trow    = row;
  tindex  = index;
  tsrow   = srow;
  tsindex = sindex;
  l_8     = (int64_t)8;
  l_4     = (int64_t)4;

  if (n > 1) {
    for (step = 1;step < n; step += step) {
      next_step = step + step;
      for (j=0;j<(n-step); j = j + next_step) {
	l1 = step;
	l2 = n-j-step;
	if (l2 > step) {
	  l2 = step;
	}
	dcrsng_mag_merge(&trow[j],&tindex[j],&trow[j+step],&tindex[j+step],
			 &tsrow[j],&tsindex[j],l1,l2);
      }
      ln = n & (next_step-1);
      if (ln <= step) {
	if (ln > 0) {
	  move_size = ln * l_8;
	  k = n-ln;
	  memcpy(&tsrow[k],&trow[k],move_size);
	  move_size = ln * l_4;
	  memcpy(&tsindex[k],&tindex[k],move_size);
	}
      }
      dtemp   = tsrow;
      tsrow   = trow;
      trow    = dtemp;
      itemp   = tsindex;
      tsindex = tindex;
      tindex  = itemp;
    } /* end for step */
    if (trow != row) {
      move_size = n * l_8;
      memcpy(row,trow,move_size);
      move_size = n * l_4;
      memcpy(index,tindex,move_size);
    }
  }
}

