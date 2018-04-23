#include "system_includes.h"
#include "dsort.h"
#include "pairwise_sum.h"
#include "dreverse_list.h"
#include "stable_add.h"
double stable_add(int n, double *v, double *scratch) {
  /*
    Return the sum of a list of n doubles, v, added in 
    stable manner. Scratch must be of length n, and is used
    for sorting the v's so as to be able to sum from 
    smallest to largest magnitude and have at most 1 
    cancellation error. Note v is reordered by this routine.
    returns the sum of elements in v.
    

    Called by: lr11_approximate_delta_concs, lr12_approximate_delta_concs,
    Calls:     dsort,pairwise_sum, dreverse_list
    
    NB might want to replace pairwise_sum with mr_sum for a still more
    stable addition.
  */
  double zero;
  double sum;
  double sumn;
  int i;
  int first_non_neg;
  int np;
  int padi;
  zero = 0.0;
  sum  = zero;
  sumn = zero;
  dsort(n,v,scratch);
  if (v[0] >= 0.0) {
    /*
      All non-negative numbers: no cancellation.
    */
    sum = pairwise_sum(n,v,scratch);
    /*
    for (i=0;i<n;i++) {
      sum += v[i];
    }
    */
  } else if (v[n-1] < 0) {
    /*
      All negative numbers: no cancellation
      We want to add from smallest magnitude to largest.
    */
    /*
    for (i=n-1;i>=0;i--) {
      sum += v[i];
    }
    */
    dreverse_list(n,v);
    sum = pairwise_sum(n,v,scratch);
  } else {
    first_non_neg = 1;
    for (i=1;i<n;i++) {
      if (v[i] >= zero) break;
      first_non_neg++;
    }
    /*
    for (i=first_non_neg;i<n;i++) {
      sum += v[i];
    }
    */
    np = n-first_non_neg;
    sum = pairwise_sum(np,&v[first_non_neg],scratch);
    /*
    for (i=first_non_neg-1;i>=0;i--) {
      sumn += v[i];
    }
    */
    dreverse_list(first_non_neg,v);
    sumn = pairwise_sum(first_non_neg,v,scratch);
    /*
      One cancellation error can occur here:
    */
    sum += sumn;
  }
  return(sum);
}

