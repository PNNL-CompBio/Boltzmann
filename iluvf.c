#include "boltzmann_structs.h"
#include "boltzmann_cvodes_headers.h"
#include "cvodes_params_struct.h"
#include "dcrsng_mag_sort.h"
#include "isort.h"
#include "iluvf.h"
int iluvf(struct state_struct *state) {
  /*
    Compute an ILU preconditioner with value based fill.
    Number of elements allowed per factor per fill row beyond
    original number is the fill level, so a fill level of zero
    will have the same number of elements as the newtion iteraton
    matrix M.
    Called by: precondition_newton_matrix
    Calls:     dcrsng_mag_sort, isort, fprintf, fflush
  */
  struct cvodes_params_struct *cvodes_params;
  double *miter;
  double *l;             
  double *u;
  double *recip_diag_u;
  double *prec_row;
  double *frow;
  double *srow;
  double recip_di;
  double multiplier;
  int    *im;
  int    *jm;
  int    *il;
  int    *jl;
  int    *iu;
  int    *ju;
  int    *lindex;
  int    *uindex;
  int    *column_mask;
  int    *sindex;

  int lpos;
  int upos;

  int epos;
  int i;

  int j;
  int k;

  int kk;
  int lcount;

  int m;
  int ucount;

  int l1;
  int lmax;
  
  int u1;
  int umax;
  
  int nnzl;
  int nnzu;

  int ny;
  int fill;

  int success;
  int padi;

  FILE *lfp;
  FILE *efp;
  success       = 1;
  lfp           = state->lfp;
  ny            = state->nunique_molecules;
  cvodes_params = state->cvodes_params;
  fill          = cvodes_params->prec_fill;
  nnzl          = cvodes_params->nnzl;
  nnzu          = cvodes_params->nnzu;
  miter         = cvodes_params->miter_m;
  l             = cvodes_params->prec_l;
  u             = cvodes_params->prec_u;
  recip_diag_u  = cvodes_params->recip_diag_u;
  prec_row      = cvodes_params->prec_row;
  frow          = cvodes_params->frow;
  srow          = cvodes_params->srow;
  im            = cvodes_params->miter_im;
  jm            = cvodes_params->miter_jm;
  il            = cvodes_params->prec_il;
  jl            = cvodes_params->prec_jl;
  iu            = cvodes_params->prec_iu;
  ju            = cvodes_params->prec_ju;
  lindex        = cvodes_params->lindex;
  uindex        = cvodes_params->uindex;
  column_mask   = cvodes_params->column_mask;
  sindex        = cvodes_params->sindex;
  lpos = 0;
  upos = 0;
  for (i=0;i<ny;i++) {
    /*
      Form row i.
    */
    il[i]  = lpos;
    iu[i]  = upos;
    lcount = 0;
    ucount = 0;
    for (k=im[i];k<im[i+1];k++) {
      j = jm[k];
      prec_row[j] = miter[k];
      column_mask[j] = 1;
      if (j < i) {
	lindex[lcount] = j;
	lcount += 1;
      } else {
	if (j > i) {
	  uindex[ucount] = j;
	  ucount += 1;
	}
      }
    }
    lmax = lcount + fill;
    if (lmax < 0) {
      lmax = 0;
    }
    umax = ucount + fill;
    if (umax < 0) {
      umax = 0;
    }
    epos = 0;
    /*
      loop over subdiagonal elements symbolically eliminating them to
      get the ordering of elimination in the full row factoring.
    */
    while (epos < lcount) {
      k = lindex[epos];
      /*
	Loop through the elements in row k of u.
      */
      for (m = iu[k];m<iu[k+1];m++) {
	j = ju[m];
	if (j <i) {
	  if (column_mask[j] == 0) {
	    column_mask[j] = 1;
	    lindex[lcount] = j;
	    lcount += 1;
	  }
	} else {
	  if (j > i) {
	    if (column_mask[j] == 0) {
	      column_mask[j] = 1;
	      uindex[ucount] = j;
	      ucount += 1;
	    }
	  }
	}
      }
      epos += 1;
    } /* end while (epos < lcount) */
    /*
      Now we need to sort lindex[0:lcount-1] so as to eliminate elements
      in order.
    */
    if (lcount > 0) {
      isort(lcount,lindex,sindex);
      for (k=0;k<lcount;k++) {
	/*
	  Eliminate element lindex[k] from prec_row.
	*/
	kk = lindex[k];
	multiplier = 0.0 - prec_row[kk] * recip_diag_u[kk];;
	for (m = iu[kk];m<iu[kk+1];m++) {
	  j = ju[m];
	  prec_row[j] += multiplier * u[m];
	} /* end for m */
	prec_row[k] = multiplier;
      } /* end for k */
    } /* end if (lcount > 0) */
    /*
      Need to test that prec_row[i] != 0 here.
      If it is, we could just arbitrarily set it to 1 to
      allow things to continue Or could return an error?
    */
    recip_di        = 1.0/prec_row[i];
    recip_diag_u[i] = recip_di;
    /*
      now we need extract lowers and uppers, and if their count > 
      lmax / umax we need to sort by decreasing magnitude keeping
      the larges lmax / umax magnitudes and discarding the rest.
    */
    l1 = lcount;
    if (lcount > lmax) {
      /*
	Extract sub diagonal elements into a vector.
	Need an extra frow vector for this, that can also be se
	for the super diagonal elements.
      */
      for(k=0;k<lcount;k++) {
	kk = lindex[k];
	frow[k] = prec_row[kk];
      }
      /*
	Now sort frow[k],lindex[k] pairs in decreasing magintude of frow[k]
	we will need 2 scratch arrays: srow[ny], sindex[ny]
      */
      dcrsng_mag_sort(lcount,frow,lindex,srow,sindex);
      /*
	Now sort the first lmax indices in the reordered lindex
	into increasing order so as to store the l matrix
	in sorted column order.
      */
      isort(lmax,lindex,sindex);
      l1 = lmax;
    }
    /*
      Extract the first l1 entries from the sorted frow,lindex pair
      into l.
    */
    if ((lpos + l1) >= nnzl) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"iluvf too many entries in L, revisit nnzl setting in boltzmann_size_jacobian\n");
	fflush(lfp);
      }
      break;
    }
    for (k=0;k<l1;k++) {
      kk = lindex[k];
      l[lpos] = prec_row[kk];
      jl[lpos] = kk;
      lpos += 1;
    }
    /*
      Now reset the subdiagonal part of prec_row and column_mask
    */
    for (k=0;k<lcount;k++) {
      kk = lindex[k];
      column_mask[kk] = 0;
      prec_row[kk]     = 0.0;
    }
    /*
      Now extract the super diagonal part.
    */
    u1 = ucount;
    if (ucount > umax) {

      u1 = umax;
      /*
	Extract super diagonal elements into a vector.
	Need an extra frow vector for this, that can also be se
	for the super diagonal elements.
      */
      for(k=0;k<ucount;k++) {
	kk = uindex[k];
	frow[k] = prec_row[kk];
      }
      dcrsng_mag_sort(ucount,frow,uindex,srow,sindex);
    }
    isort(u1,uindex,sindex);
    /*
      Extract the first u1 entries from the sorted frow,uindex pair
      into u. 
    */
    if (upos + u1 > nnzu) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"iluvf too many entries in U, revisit nnzu setting in boltzmann_size_jacobian\n");
	fflush(lfp);
      }
      break;
    }
    for (k=0;k<u1;k++) {
      kk = uindex[k];
      u[upos] = prec_row[kk];
      ju[upos] = kk;
      upos += 1;
    }
    /*
      Now reset the super diagonal part of prec_row and column_mask
    */
    for (k=0;k<ucount;k++) {
      kk = uindex[k];
      column_mask[kk] = 0;
      prec_row[kk]     = 0.0;
    }
  } /* end for (i...) */
  il[ny] = lpos;
  iu[ny] = upos;
  return(success);
}
