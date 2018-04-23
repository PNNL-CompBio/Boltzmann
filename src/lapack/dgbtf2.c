#include "boltzmann_structs.h"
#include "blas.h"
#include "dger.h"
#include "dswap.h"
#include "dgbtf2.h"
/*
*> \brief \b DGBTF2 computes the LU factorization of a general band matrix using the unblocked version of the algorithm.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DGBTF2 + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbtf2.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbtf2.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbtf2.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
* 
*       .. Scalar Arguments ..
*       INTEGER            INFO, KL, KU, LDAB, M, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   AB( LDAB, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGBTF2 computes an LU factorization of a real m-by-n band matrix A
*> using partial pivoting with row interchanges.
*>
*> This is the unblocked version of the algorithm, calling Level 2 BLAS.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
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
*> \param[in,out] AB
*> \verbatim
*>          AB is DOUBLE PRECISION array, dimension (LDAB,N)
*>          On entry, the matrix A in band storage, in rows KL+1 to
*>          2*KL+KU+1; rows 1 to KL of the array need not be set.
*>          The j-th column of A is stored in the j-th column of the
*>          array AB as follows:
*>          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
*>
*>          On exit, details of the factorization: U is stored as an
*>          upper triangular band matrix with KL+KU superdiagonals in
*>          rows 1 to KL+KU+1, and the multipliers used during the
*>          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
*>          See below for further details.
*> \endverbatim
*>
*> \param[in] LDAB
*> \verbatim
*>          LDAB is INTEGER
*>          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (min(M,N))
*>          The pivot indices; for 1 <= i <= min(M,N), row i of the
*>          matrix was interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*>          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
*>               has been completed, but the factor U is exactly
*>               singular, and division by zero will occur if it is used
*>               to solve a system of equations.
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
*> \date September 2012
*
*> \ingroup doubleGBcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The band storage scheme is illustrated by the following example, when
*>  M = N = 6, KL = 2, KU = 1:
*>
*>  On entry:                       On exit:
*>
*>      *    *    *    +    +    +       *    *    *   u14  u25  u36
*>      *    *    +    +    +    +       *    *   u13  u24  u35  u46
*>      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
*>     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
*>     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
*>     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
*>
*>  Array elements marked * are not used by the routine; elements marked
*>  + need not be set on entry, but are required by the routine to store
*>  elements of U, because of fill-in resulting from the row
*>  interchanges.
*> \endverbatim
*>
*  =====================================================================

      SUBROUTINE DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
*/
void dgbtf2_(int *mp ,int *np, int *klp, int *kup, double *ab_p,
	     int *ldabp, int *ipiv_p, int *infop) {
  /*
*
*  -- LAPACK computational routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   AB( LDAB, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
*/
  double *ab;
  double one;
  double zero;
  double m_one;
  double multiplier;
  int *ipiv;

  int m;
  int n;

  int kl;
  int ku;

  int ldab;
  int info;

  int i;
  int j;

  int jp;
  int ju;

  int km;
  int kv;

  int min_kv_n;
  int min_m_n;

  int kv_ldab;
  int j_ldab;
  
  int kmp1;
  int i_one;

  int ju_cand;
  int jumjp1;

  int jumj;
  int jpkv_ldab;

  int ldabm1;
  int padi;

  one = 1.0;
  zero = 0.0;
  m_one = zero-one;
  i_one = 1;
  m    = *mp;
  n    = *np;
  kl   = *klp;
  ku   = *kup;
  ldab  = *ldabp;
  ldabm1 = ldab - 1;
  /* 
    Caution address arithmetic here in computing ab and ipiv addresses
    Used to avoid having to subtract ldab from the second dimension
    and 1 from the first dimension in two dimensional references to ab,
    and to avoid subtraction of 1 from ipiv indices, arrays in c start
    at 0
  */
  ab   = ab_p - (1 + ldab);
  ipiv = ipiv_p - 1;
  /*  

      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, JP, JU, KM, KV
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      EXTERNAL           IDAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     KV is the number of superdiagonals in the factor U, allowing for
*     fill-in.
*
      KV = KU + KL
  */
  kv = ku + kl;
  /*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGBTF2', -INFO )
         RETURN
      END IF
  */
  info = 0;
  if (m < 0) {
    info = -1;
  } else if (n < 0) {
    info = -2;
  } else if (kl < 0) {
    info = -3;
  } else if (ku < 0) {
    info = -4;
  } else if (ldab < (kl + kv + 1)) {
    info = -6;
  }
  if (info != 0) {
    fprintf(stderr,"dgbtf2: Error info = %d\n",info);
    fflush(stderr);
  } else {
    /*    
      Quick return if possible
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
    */
    if ((m > 0) && (n > 0)) {
      /*
*     Gaussian elimination with partial pivoting
*
*     Set fill-in elements in columns KU+2 to KV to zero.
*
      DO 20 J = KU + 2, MIN( KV, N )
         DO 10 I = KV - J + 2, KL
            AB( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
      */
      min_kv_n = kv;
      if (n < kv) {
	min_kv_n = n;
      }
      j_ldab = ldab * (ku + 2);
      kv_ldab = ldab * kv;
      for (j= ku + 2; j <= min_kv_n;j++) {
	for (i=kv-j+2;i<=kl;i++) {
	  ab[j_ldab + i] = zero;
	}
	j_ldab += ldab;
      }
      /*					
*
*     JU is the index of the last column affected by the current stage
*     of the factorization.
*
      JU = 1
      */
      ju = 1;
      /*
*
      DO 40 J = 1, MIN( M, N )
*
*        Set fill-in elements in column J+KV to zero.
*
         IF( J+KV.LE.N ) THEN
            DO 30 I = 1, KL
               AB( I, J+KV ) = ZERO
   30       CONTINUE
         END IF
      */
      min_m_n = m;
      if (n < m) {
	min_m_n = n;
      }
      j_ldab = ldab;
      jpkv_ldab = ldab + kv_ldab;
      for (j=1;j<= min_m_n;j++) {
	if (j+kv <= n) {
	  for ( i=1;i<=kl;i++) {
	    ab[jpkv_ldab + i] = zero;
	  }
	}
	/*
*        Find pivot and test for singularity. KM is the number of
*        subdiagonal elements in the current column.
         KM = MIN( KL, M-J )
         JP = IDAMAX( KM+1, AB( KV+1, J ), 1 )
         IPIV( J ) = JP + J - 1
        */
	km = kl;
	if ((m-j) < km) {
	  km = m-j;
	}
	kmp1 = km + 1;
	jp = idamax_(&kmp1, &ab[j_ldab + kv+1], &i_one);
	ipiv[j] = jp + j-1;
	/*	
         IF( AB( KV+JP, J ).NE.ZERO ) THEN
            JU = MAX( JU, MIN( J+KU+JP-1, N ) )
*
*           Apply interchange to columns J to JU.
*
            IF( JP.NE.1 )
     $         CALL DSWAP( JU-J+1, AB( KV+JP, J ), LDAB-1,
     $                     AB( KV+1, J ), LDAB-1 )
*
            IF( KM.GT.0 ) THEN
*
*              Compute multipliers.
*
               CALL DSCAL( KM, ONE / AB( KV+1, J ), AB( KV+2, J ), 1 )
*
*              Update trailing submatrix within the band.
*
               IF( JU.GT.J )
     $            CALL DGER( KM, JU-J, -ONE, AB( KV+2, J ), 1,
     $                       AB( KV, J+1 ), LDAB-1, AB( KV+1, J+1 ),
     $                       LDAB-1 )
            END IF
	*/
	if (ab[j_ldab + kv + jp] != 0) {
	  ju_cand = j + ku + jp-1;
	  if (ju_cand > n) {
	    ju_cand = n;
	  } 
	  if (ju_cand > ju) {
	    ju = ju_cand;
	  }
	  jumj = ju - j;
	  if (jp != 1) {
	    jumjp1 = jumj + 1;
	    dswap_(&jumjp1,&ab[j_ldab + kv +jp ],&ldabm1,
		   &ab[j_ldab+kv+1],&ldabm1);
	  }
	  if (km > 0) {
	    multiplier = one/ab[j_ldab+kv+1];
	    dscal_(&km, &multiplier, &ab[j_ldab+kv+2], &i_one);
	    if (ju > j) {
	      dger_(&km,&jumj,&m_one,&ab[j_ldab+kv+2],&i_one,
		    &ab[j_ldab+ldab+kv],&ldabm1,&ab[j_ldab+ldab+kv+1],
		    &ldabm1);
	    }
	  }
	} else {
	  /*
         ELSE
           If pivot is zero, set INFO to the index of the pivot
           unless a zero pivot has already been found.
            IF( INFO.EQ.0 )
     $         INFO = J
         END IF

	  */
	  if (info == 0) {
	    info = j;
	  }
	}
	j_ldab += ldab;
	jpkv_ldab += ldab;
	/*
   40 CONTINUE
   
      RETURN
	*/
      } /* end for (j...) */
      /*
*
*     End of DGBTF2
*
      END
      */
    } /* end if (im > 0 and in > 0) */
  } /* end else info was 0 */
}
