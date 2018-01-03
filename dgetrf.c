#include "system_includes.h"
#include "dgetrf2.h"
#include "dlaswp.h"
#include "dgemm.h"
#include "dtrsm.h"
#include "dgetrf.h"
void dgetrf_(int *m_p, 
	    int *n_p, 
	    double *a, 
	    int *lda_p, 
	    int *ipiv, 
	    int *info_p) {
  /*
    Adapted from dgetrf.f from lapack source.
  */
  double one;
  double mone;
  double zero;
  double temp;
  double sfmin;
  double *dbl_ptr;
  int64_t sfmin_hex;
  int64_t ajj_pos;
  int64_t jb_lda;
  int64_t j_lda;
  int64_t j_p_jb_lda;
  int64_t a21_start;
  int64_t a22_start;
  int m;
  int n;

  int lda;
  int info;

  int iinfo;
  int i;

  int j;
  int jb;
 
  int nb;
  int inc1;

  int izero;
  int min_m_n;

  int mmj;
  int padi;

  int min_m_j_p_jb;
  int j_p_jb_m1;

  int n_m_j_m_jb;
  int m_m_j_m_jb;

  
  char l_char;
  char n_char;
  char u_char;
  char padc;
  m   = *m_p;
  n   = *n_p;
  lda = *lda_p;
  sfmin_hex       = 0x0004000000000001L;
  dbl_ptr = (double*)&sfmin_hex;
  sfmin   = *dbl_ptr;
  one     = 1.0;
  zero    = 0.0;
  mone    = -1.0;
  inc1    = 1;
  l_char  = 'L';
  n_char  = 'N';
  u_char  = 'U';
  nb = 64;
  info = 0;
  if (m < 0) {
    info = -1;
  } else {
    if (n < 0) {
      info = -2;
    } else {
      if ((lda < 1) || (lda < m)) {
	info = -4;
      }
    }
  }
  if (info != 0) {
    fprintf(stderr,"dgetrf: error: info = %d\n",info);
    fflush(stderr);
  } else {
    if ((m > 0) && (n > 0)) {
      min_m_n = (n < m) ? n : m;
      if ((nb <= 1) || (nb >=  min_m_n)) {
	/*
	  Use point algorithm.
	*/
	dgetrf2_(&m, &n, a, &lda, ipiv, &info);
      } else {
	/*
	  Use block algorithm.
	*/
	j_lda = 0;
	jb_lda = nb * lda;
	for (j = 0;j < min_m_n;j += nb) {
	  jb = ((min_m_n -j) < nb) ? (min_m_n - j) : nb;
	  jb_lda = jb * lda;
	  j_p_jb_lda = j_lda + jb_lda;
	  /*
	    Factor diagonal and subdiagonal blocks and test for exact
	    singularity.
	  */
	  mmj = m-j;
	  ajj_pos = j_lda + j;
          dgetrf2_(&mmj, &jb, &a[ajj_pos], &lda, &ipiv[j], &iinfo);
	  /*
            Adjust INFO and the pivot indices.
	  */
	  if ((info == 0) && (iinfo > 0) ) {
	    info = iinfo + j;
	  }
	  min_m_j_p_jb = (m < (j+jb)) ? m : (j+jb) ;
	  for (i=j;i<min_m_j_p_jb;i++) {
	    /*
	      Here we don't mess with the start from one
	      fortranism for ipiv as it is on both sides of the
	      assignment.
	    */
	    ipiv[i] = j + ipiv[i];
	  }
	  j_p_jb_m1 = j + jb - 1;
	  /*
            Apply interchanges to columns 0:j-1.
	  */
          dlaswp_( &j, a, &lda, &j, &j_p_jb_m1, ipiv, &inc1 );
	  
	  n_m_j_m_jb = n - j - jb;
	  m_m_j_m_jb = m - j - jb;
	  if (n_m_j_m_jb > 0) {
	    /*
	      Apply interchanges to columns J+JB:N-1.
	    */
	    dlaswp_(&n_m_j_m_jb, &a[j_p_jb_lda], &lda, &j, &j_p_jb_m1,
		   ipiv, &inc1);
	    /*
	      Compute block row of U.
	    */
	    dtrsm_(&l_char,&l_char, &n_char, &u_char,
		  &jb, &n_m_j_m_jb, &one, &a[ajj_pos], &lda,
		  &a[j+j_p_jb_lda], &lda);

	    if (m_m_j_m_jb > 0) {
	      /*
                Update trailing submatrix.
	      */
	      dgemm_(&n_char, &n_char, &m_m_j_m_jb,
		    &n_m_j_m_jb, &jb,  &mone,
		    &a[j_lda + j + jb],&lda,
		    &a[j_p_jb_lda + j],&lda,&one,
		    &a[j_p_jb_lda + j + jb],&lda);
	    } /* end if more rows left */
	  } /* end if more columns left. */
	  j_lda = j_p_jb_lda;
	} /* end for j */
      } /* end block algorithm */
    } /* end nontrivial system */
  } /* end else arguments ok */
  *info_p = info;
} /* end dgetrf */
