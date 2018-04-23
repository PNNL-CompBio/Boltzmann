/* dgbtrs.c
*> \brief \b DGBTRS
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DGBTRS + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbtrs.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbtrs.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbtrs.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB,
*                          INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          TRANS
*       INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGBTRS solves a system of linear equations
*>    A * X = B  or  A**T * X = B
*> with a general band matrix A using the LU factorization computed
*> by DGBTRF.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          Specifies the form of the system of equations.
*>          = 'N':  A * X = B  (No transpose)
*>          = 'T':  A**T* X = B  (Transpose)
*>          = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] KL
*> \verbatim
*>          KL is INTEGER
*>          The number of subdiagonals within the band of A.  KL >= 0.
*> \endverbatim
*>
*> \param[in] KU
*> \verbatim
*>          KU is INTEGER
*>          The number of superdiagonals within the band of A.  KU >= 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of right hand sides, i.e., the number of columns
*>          of the matrix B.  NRHS >= 0.
*> \endverbatim
*>
*> \param[in] AB
*> \verbatim
*>          AB is DOUBLE PRECISION array, dimension (LDAB,N)
*>          Details of the LU factorization of the band matrix A, as
*>          computed by DGBTRF.  U is stored as an upper triangular band
*>          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
*>          the multipliers used during the factorization are stored in
*>          rows KL+KU+2 to 2*KL+KU+1.
*> \endverbatim
*>
*> \param[in] LDAB
*> \verbatim
*>          LDAB is INTEGER
*>          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          The pivot indices; for 1 <= i <= N, row i of the matrix was
*>          interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
*>          On entry, the right hand side matrix B.
*>          On exit, the solution matrix X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup doubleGBcomputational
*
*  =====================================================================
*/
/*
      SUBROUTINE DGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB,
     $                   INFO )
*/
#include "system_includes.h"
#include "blas.h"
#include "dger.h"
#include "dswap.h"
#include "dtbsv.h"
#include "lsame.h"
#include "dgbtrs.h"
void dgbtrs_(char *transp, int *np, int *klp, int *kup, int *nrhsp,
	     double *abp, int *ldabp, int *ipivp, double *bp, int *ldbp, 
	     int *infop) {
/*
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*/
  double *ab;
  double *b;
  double one;
  double m_one;
  int  *ipiv;

  int  info;
  int  n;

  int  kl;
  int  ku;

  int  nrhs;
  int  ldab;

  int  ldb;
  int  lnoti;

  int  notran;
  int  tran;

  int  conj;
  int  valid_arg1;

  int  i;
  int  j;

  int  kd;
  int  l;

  int  lm;
  int  i_one;
  
  int j_ldab;
  int kl_p_ku;

  int i_ldb;
  int padi;

  char trans;
  char uchar;
  char nchar;
  char tchar;
  char cchar;
  char cpad1;
  char cpad2;
  char cpad3;
/*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LNOTI, NOTRAN
      INTEGER            I, J, KD, L, LM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DGER, DSWAP, DTBSV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*/
  one = 1.0;
  m_one = -1.0;
  i_one = 1;
  n    	= *np;
  kl   	= *klp;
  ku   	= *kup;
  nrhs 	= *nrhsp;
  ldab 	= *ldabp;
  ldb  	= *ldbp;
  trans = transp[0];
  /*
    Caution address arthmetic in the next three statments to 
    account for difference in fortran and c array indices.
  */
  ab   	= abp - 1 - ldab;
  b    	= bp - 1 - ldb;
  ipiv 	= ipivp - 1;
  /*  
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
  */
  uchar  = 'U';
  nchar  = 'N';
  tchar  = 'T';
  cchar  = 'C';
  kl_p_ku = kl + ku;
  notran = lsame_(&trans,&nchar,1,1);
  tran   = lsame_(&trans,&tchar,1,1);
  conj   = lsame_(&trans,&cchar,1,1);
  valid_arg1 = notran || tran || conj;
  info = 0;
  if (! valid_arg1) {
    info = -1;
  } else if (n < 0) {
    info = -2;
  } else if (kl < 0) {
    info = -3;
  } else if (ku < 0) {
    info = -4;
  } else if (nrhs < 0) {
    info = -5;
  } else if (ldab < ((2*kl) + ku + 1)) {
    info = -7;
  } else if ((ldb < n) || (ldb < 1)) {
    info = -10;
  }
  /*
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.( 2*KL+KU+1 ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGBTRS', -INFO )
         RETURN
      END IF
  */
  if (info == 0) {
    /*    
*     Quick return if possible

      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
*/
    if ((n > 0) && (nrhs > 0) ) {
      /*
      KD = KU + KL + 1
      LNOTI = KL.GT.0
      */
      kd = ku + kl + 1;
      lnoti = kl > 0;
      /*
      IF( NOTRAN ) THEN
         IF( LNOTI ) THEN
*
*        Solve  A*X = B.
*
*        Solve L*X = B, overwriting B with X.
*
*        L is represented as a product of permutations and unit lower
*        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
*        where each transformation L(i) is a rank-one modification of
*        the identity matrix.
*
*/
      if (notran) {
	if (lnoti) {
	  /*
            DO 10 J = 1, N - 1
               LM = MIN( KL, N-J )
               L = IPIV( J )
               IF( L.NE.J )
     $            CALL DSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
               CALL DGER( LM, NRHS, -ONE, AB( KD+1, J ), 1, B( J, 1 ),
     $                    LDB, B( J+1, 1 ), LDB )
   10       CONTINUE
         END IF
	  */
	  j_ldab = ldab;
	  for (j=1;j<=n-1;j++) {
	    lm = n-j;
	    if (kl < lm) {
	      lm = kl;
	    }
	    l = ipiv[j];
	    if (l != j) {
	      dswap_(&nrhs,&b[ldb + l],&ldb,&b[ldb + j],&ldb);
	    }
	    dger_(&lm,&nrhs,&m_one,&ab[j_ldab + kd + 1],&i_one,
		  &b[j+ldb],&ldb,&b[j+1+ldb],&ldb);
	    j_ldab += ldab;
	  } /* end for j */
	} /* end if lnoti */
	      
	/*
         DO 20 I = 1, NRHS
*
*           Solve U*X = B, overwriting B with X.
*
            CALL DTBSV( 'Upper', 'No transpose', 'Non-unit', N, KL+KU,
     $                  AB, LDAB, B( 1, I ), 1 )
   20    CONTINUE
	*/
	i_ldb = ldb;
	for (i=1;i<=nrhs;i++) {
	  dtbsv_(&uchar,&nchar,&nchar,&n,&kl_p_ku,abp,&ldab,&b[1+i_ldb],&i_one);
	  i_ldb += ldb;
	}
      } else {
	/*	    
      ELSE
*
*        Solve A**T*X = B.
*
         DO 30 I = 1, NRHS
*
*           Solve U**T*X = B, overwriting B with X.
*
            CALL DTBSV( 'Upper', 'Transpose', 'Non-unit', N, KL+KU, AB,
     $                  LDAB, B( 1, I ), 1 )
   30    CONTINUE
	*/
	i_ldb = ldb;
	for (i=1;i<=nrhs;i++) {
	  dtbsv_(&uchar,&tchar,&nchar,&n,&kl_p_ku,abp,&ldab,&b[1+i_ldb],&i_one);
	  i_ldb += ldb;
	}
	/*
          Solve L**T*X = B, overwriting B with X.

         IF( LNOTI ) THEN
            DO 40 J = N - 1, 1, -1
               LM = MIN( KL, N-J )
               CALL DGEMV( 'Transpose', LM, NRHS, -ONE, B( J+1, 1 ),
     $                     LDB, AB( KD+1, J ), 1, ONE, B( J, 1 ), LDB )
               L = IPIV( J )
               IF( L.NE.J )
     $            CALL DSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
   40       CONTINUE
         END IF
	*/
	if (lnoti) {
	  j_ldab= (n-1) * ldab;
	  for (j=n-1;j>=1;j--) {
	    lm = n-j;
	    if (kl < lm) {
	      lm = kl;
	    }
	    dgemv_(&tchar,&lm,&nrhs,&m_one,&b[j+1+ldb],&ldb,
		   &ab[j_ldab + kd + 1],&i_one,&one,&b[j+ldb],&ldb,1);
	    l = ipiv[j];
	    if (l != j) {
	      dswap_(&nrhs,&b[l+ldb],&ldb,&b[j+ldb],&ldb);
	    }
	    j_ldab -= ldab;
	  }
	}
	/*
      END IF
      RETURN

      End of DGBTRS
      END
	*/
      } /* end else solve transpose problem. */
    } /* end not quick returnt */
  } /* end if valid arguments */
}
