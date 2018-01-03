#include "system_includes.h"
#include "lsame.h"
#include "dlaswp.h"
#include "dtrsm.h"
#include "dgetrs.h"
void dgetrs(char *trans_p, 
	    int *n_p, 
	    int *nrhs_p, 
	    double *a,
	    int *lda_p,
	    int *ipiv,
	    double *b,
	    int *ldb_p,
	    int *info_p) {
  /*
    Adapted from dgetrs.f from netlib.
    Solves system of linear equations

       a x = b or a^t x = b with 
       n by n matrix a using the lu factorization computed by dgetrf.

    Calls: lsame, dlaswp, dtrsm, fprintf, fflush
    
  */
  /*
    Test the input parmaeters.
  */
  double one;
  double zero;
  int  n;
  int  nrhs;

  int  lda;
  int  ldb;

  int  info;
  int  notran;

  int  izero;
  int  inc1;

  int  nm1;
  int  incm1;
    
  char trans;
  char n_char;
  char t_char;
  char c_char;

  char l_char;
  char r_char;
  char u_char;

  n = *n_p;
  nrhs = *nrhs_p;
  lda  = *lda_p;
  ldb  = *ldb_p;
  trans = *trans_p;
  n_char = 'N';
  t_char = 'T';
  c_char = 'C';
  l_char = 'L';
  r_char = 'R';
  u_char = 'U';
  nm1    = n-1;
  one  = 1.0;
  zero = 0.0;
  info = 0;
  izero = 0;
  inc1  = 1;
  incm1 = -1;
  notran = lsame(&trans,&n_char);
  if ((notran == 0) && ((lsame(&trans,&t_char) + lsame(&trans,&c_char)) == 0)) {
    info = -1;
  } else {
    if (n < 0) {
      info = -2;
    } else {
      if (nrhs < 0) {
	info = -3;
      } else {
	if ((lda < 1) || (lda < n)) {
	  info = -5;
	} else {
	  if ((ldb < 1) || (ldb < n)) {
	    info = -8;
	  }
	}
      }
    }
  }
  if (info != 0) {
    fprintf (stderr,"dgetrs: Error in arguments info = %d\n",info);
    fflush(stderr);
  } else {
    if ((n > 0) && (nrhs > 0)) {
      if (notran) {
	/*
	  Solve  a * x = b.
	*/
	/*
	  Apply row interchanges to the right hand sidse.
	*/
	dlaswp(&nrhs, b, &ldb, &izero, &nm1, ipiv, &inc1);
	/*
	  Solve L x = b overwriting B with X
	*/
	dtrsm(&l_char,&l_char,&n_char,&u_char,&n,&nrhs,
	      &one, a, &lda, b, &ldb);
	/*
	  Solve U x = b overwriting b with x.
	*/
	dtrsm(&l_char,&u_char,&n_char,&n_char,&n,&nrhs,
	      &one,a, &lda, b, &ldb);
      }
    } else {
      /*
	solve a^t * x = b
      */
      /* 
	 solve u^t *x = b overwriting b with x
      */
      dtrsm(&l_char,&u_char,&t_char,&n_char,&n,&nrhs,
	    &one,a,&lda,b,&ldb);
      /*
	solve l^t * x = b overwriting b with x
      */
      dtrsm(&l_char,&l_char,&t_char,&u_char,&n,&nrhs,
	    &one, a,&lda,b,&ldb);
      /*
	Apply row interchanges to solution vectors.
      */
      dlaswp(&nrhs, b, &ldb, &izero, &nm1, ipiv, &incm1);
    } /* end else notran == 0 */
  } /* end else info == 0 */
  *info_p = info;
}
