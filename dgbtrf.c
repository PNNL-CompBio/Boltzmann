#include "system_includes.h"
#include "blas.h"
#include "dgbtf2.h"
#include "dlaswp.h"
#include "dger.h"
#include "dgemm.h"
#include "dtrsm.h"
#include "dswap.h"
#include "dgbtrf.h"

void dgbtrf_(int *mp, int *np, int *klp, int *kup, double *ab_p, int *ldabp,
	     int *ipiv_p, int *infop) {
/*
*
*       SUBROUTINE DGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
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
*> DGBTRF computes an LU factorization of a real m-by-n band matrix A
*> using partial pivoting with row interchanges.
*>
*> This is the blocked version of the algorithm, calling Level 3 BLAS.
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
v*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
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
*>  elements of U because of fill-in resulting from the row interchanges.
*> \endverbatim
*>
*  =====================================================================
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      INTEGER            NBMAX, LDWORK
      PARAMETER          ( NBMAX = 64, LDWORK = NBMAX+1 )

*     .. Local Scalars ..
      INTEGER            I, I2, I3, II, IP, J, J2, J3, JB, JJ, JM, JP,
     $                   JU, K2, KM, KV, NB, NW
      DOUBLE PRECISION   TEMP

*     .. Local Arrays ..
      DOUBLE PRECISION   WORK13( LDWORK, NBMAX ),
     $                   WORK31( LDWORK, NBMAX )

*     .. External Functions ..
      INTEGER            IDAMAX, ILAENV
      EXTERNAL           IDAMAX, ILAENV

     .. External Subroutines ..
      EXTERNAL           DCOPY, DGBTF2, DGEMM, DGER, DLASWP, DSCAL,
     $                   DSWAP, DTRSM, XERBLA

*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
     ..
*/
  double *ab;
  double *work13;
  double *work31;
  double one;
  double m_one;
  double zero;
  double work13_p[4160];
  double work31_p[4160];
  double temp;
  double multiplier;

  int *ipiv;
  int m;
  int n;

  int kl;
  int ku;

  int ldab;
  int info;

  int nbmax;
  int ldwork;

  int i;
  int i2;

  int i3;
  int ip;

  int j;
  int j2;

  int j3;
  int jb;

  int jj;
  int jm;

  int jp;
  int ju;

  int k2;
  int km;

  int kv;
  int nb;

  int min_kv_n;
  int min_m_n;

  int kmp1;  
  int ione;

  int ju_cand;
  int ldabm1;

  int jjmj;
  int jmjjpjb;

  int j_ldab;
  int jj_ldab;

  int jm_m_jj;
  int jb_ldab;

  int jpjb_ldab;
  int jj_ldwork;

  int jjpjpkv_ldab;
  int kv_ldab;

  int jpkv_ldab;
  int j_ldwork;

  int row_pos;
  int ii;

  int nw;
  int padi;

  char l_char;
  char n_char;
  char u_char;
  char padc;


  one = 1.0;
  zero = 0.0;
  m_one = zero - one;
  ione = 1;
  nbmax = 64;
  ldwork = 65;
  l_char  = 'L';
  n_char  = 'N';
  u_char  = 'U';

  m    = *mp;
  n    = *np;
  kl   = *klp;
  ku   = *kup;
  ldab = *ldabp;
  /* 
    Caution address arithmetic here in computing ab and ipiv addresses
    Used to avoid having to subtract ldab from the second dimension
    and 1 from the first dimension in two dimensional references to ab,
    and to avoid subtraction of 1 from ipiv indices, arrays in c start
    at 0, similarly for work13 and work31
  */
  ab   = ab_p - (ldab + 1);
  ipiv = ipiv_p - 1;
  work13 = work13_p - (ldwork + 1);
  work31 = work31_p - (ldwork + 1);

  /*
    KV is the number of superdiagonals in the factor U, allowing for
    fill-in
  */
  kv = ku + kl;
  kv_ldab = kv * ldab;
  /*
    Test the input parameters.
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
  } else if (ldab < (kl + ku + 1)) {
    info = -6;
  }
  if (info != 0) {
    fprintf(stderr,"dgbtrf: Error, info = %d\n",info);
    fflush(stderr);
  } 
  /*
    Quick return if possible
  */
  if ((info == 0) && (m > 0) && (n > 0)) {
    /*
      Determine the block size for this environment
      
      NB = ILAENV( 1, 'DGBTRF', ' ', M, N, KL, KU )
    */
    if (ku <= 64) {
      nb = 1;
    } else {
      nb = 32;
    }
    /*
     *     The block size must not exceed the limit set by the size of the
     *     local arrays WORK13 and WORK31.
     *
      NB = MIN( NB, NBMAX )
    */
    if ((nb <= 1) || (nb > kl)) {
      /*
        Use unblocked code
      */
      dgbtf2_(mp,np,klp,kup, ab,ldabp,ipiv,infop);
    } else {
      /*
	Use blocked code

        Zero the superdiagonal elements of the work array WORK13
      */
      /*
         DO 20 J = 1, NB
            DO 10 I = 1, J - 1
               WORK13( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      */
      j_ldwork = ldwork;
      for (j=1;j<=nb;j++) {
	for (i=1;i<=j-1;i++) {
	  work13[j_ldwork +i] = zero;
	  /*
	  work13[i,j] = zero;
	  */
	}
	j_ldwork += ldwork;
      }
      /*
        Zero the subdiagonal elements of the work array WORK31
      */
      /*
         DO 40 J = 1, NB
            DO 30 I = J + 1, NB
               WORK31( I, J ) = ZERO
   30       CONTINUE
   40    CONTINUE
      */
      j_ldwork = ldwork;
      for (j=1;j<=nb;j++) {
	for (i=j+1;i<=nb;i++) {
	  work31[j_ldwork+i] = 0;
	}
	j_ldwork += ldwork;
      }
      /*
        Gaussian elimination with partial pivoting

        Set fill-in elements in columns KU+2 to KV to zero
      */
      min_kv_n = kv;
      if (n < kv) {
	min_kv_n = n;
      }
      /*
         DO 60 J = KU + 2, MIN( KV, N )
            DO 50 I = KV - J + 2, KL
               AB( I, J ) = ZERO
   50       CONTINUE
   60    CONTINUE
      */
      j_ldab = (ku+2) * ldab;
      for (j=ku+2;j<=min_kv_n;j++) {
	for (i = kv - j + 2;i<=kl;i++) {
	  ab[j_ldab + i] = zero;
	}
	j_ldab += ldab;
      }
      /*
        JU is the index of the last column affected by the current
        stage of the factorization
	JU = 1;
      */
      ju = 1;
      /*
         DO 180 J = 1, MIN( M, N ), NB
            JB = MIN( NB, MIN( M, N )-J+1 )
      */
      min_m_n = m;
      if (n < min_m_n) {
	min_m_n = n;
      }
      j_ldab = ldab;
      for (j=1;j<=min_m_n;j += nb) {
	jb = nb;
	if ((min_m_n - j+1) < jb) {
	  jb = min_m_n - j + 1;
	}
	jb_ldab = jb * ldab;
	jpjb_ldab = jb_ldab + j_ldab;
	/*
	  The active part of the matrix is partitioned
	  
              A11   A12   A13
              A21   A22   A23
              A31   A32   A33

           Here A11, A21 and A31 denote the current block of JB columns
           which is about to be factorized. The number of rows in the
           partitioning are JB, I2, I3 respectively, and the numbers
           of columns are JB, J2, J3. The superdiagonal elements of A13
           and the subdiagonal elements of A31 lie outside the band.
            I2 = MIN( KL-JB, M-J-JB+1 )
            I3 = MIN( JB, M-J-KL+1 )
	*/
	i2 = kl-jb;
	if ((m-j-jb+1) < i2) {
	  i2 = m-j-jb+1;
	}
	i3 = jb;
	if ((m-j-kl+1) < i3) {
	  i3 = m-j-kl+1;
	}
	/*

          J2 and J3 are computed after JU has been updated.

          Factorize the current block of JB columns
          DO 80 JJ = J, J + JB - 1
	*/
	jj_ldab = j_ldab;
	for (jj = j;jj<=j+jb-1;jj++) {
	  jjmj = jj - j;
	  /*
	    Set fill-in elements in column JJ+KV to zero
               IF( JJ+KV.LE.N ) THEN
                  DO 70 I = 1, KL
                     AB( I, JJ+KV ) = ZERO
   70             CONTINUE
               END IF
	  */
	  if (jj + kv <= n) {
	    row_pos = jj_ldab+kv_ldab;
	    for (i=1;i<=kl;i++) {
	      ab[row_pos+i] = zero;
	    }
	  }
	  /*
            Find pivot and test for singularity. KM is the number of
            subdiagonal elements in the current column.
	    KM = MIN( KL, M-JJ )
	  */
	  km = kl;
	  if (m-jj < km) {
	    km = m - jj;
	  }
	  /*
	    JP = IDAMAX( KM+1, AB( KV+1, JJ ), 1 )
	    IPIV( JJ ) = JP + JJ - J
	  */

	  kmp1 = km + 1;
	  jp = idamax_(&kmp1,&ab[jj_ldab + kv+1],&ione);
	  ipiv[jj] = jp + jjmj;
	  /*
               IF( AB( KV+JP, JJ ).NE.ZERO ) THEN
                  JU = MAX( JU, MIN( JJ+KU+JP-1, N ) )
                  IF( JP.NE.1 ) THEN
	  */
	  if (ab[jj_ldab+kv+jp] != zero) {
	    ju_cand = (jj + ku + jp-1);
	    if (n < ju_cand) {
	      ju_cand = n;
	    }
	    if (ju_cand > ju) {
	      ju = ju_cand;
	    }
	    if (jp != 1) {
	      /*

		Apply interchange to columns J to J+JB-1
                     IF( JP+JJ-1.LT.J+KL ) THEN
                        CALL DSWAP( JB, AB( KV+1+JJ-J, J ), LDAB-1,
		                    AB( KV+JP+JJ-J, J ), LDAB-1 )
	      */
	      if ((jp + jj - 1)  < (j + kl)) {
		ldabm1 = ldab - 1;
		dswap_(&jb,&ab[j_ldab + kv + 1 + jj - j],&ldabm1,
		       &ab[j_ldab + kv + jp + jj-j],&ldabm1);
	      } else {
		/*
                     ELSE

                       The interchange affects columns J to JJ-1 of A31
                       which are stored in the work array WORK31
                        CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,
     $                              WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                        CALL DSWAP( J+JB-JJ, AB( KV+1, JJ ), LDAB-1,
     $                              AB( KV+JP, JJ ), LDAB-1 )
                     END IF
                  END IF

		*/
		jjmj = jj - j;
		jmjjpjb = jb - jjmj;
		dswap_(&jjmj,&ab[j_ldab + kv +1 + jjmj],&ldabm1,
		       &work31[ldwork + jp+jjmj-kl],&ldwork);

		dswap_(&jmjjpjb,&ab[kv + 1 + jj_ldab], &ldabm1,
		       &ab[kv+jp + jj_ldab],&ldabm1);
	      }
	    } /* end if (jp != 0) */
	    /*
	      Compute multipliers
	      CALL DSCAL( KM, ONE / AB( KV+1, JJ ), AB( KV+2, JJ ),
     $                        1 )
	    */
	    multiplier = one / ab[jj_ldab + kv + 1];
	    dscal_(&km,&multiplier,&ab[jj_ldab+kv+2],&ione);
	    /*
	      Update trailing submatrix within the band and within
              the current block. JM is the index of the last column
              which needs to be updated.
                  JM = MIN( JU, J+JB-1 )
                  IF( JM.GT.JJ )
     $               CALL DGER( KM, JM-JJ, -ONE, AB( KV+2, JJ ), 1,
     $                          AB( KV, JJ+1 ), LDAB-1,
     $                          AB( KV+1, JJ+1 ), LDAB-1 )
	    */
	    jm = ju;
	    if (j + jb - 1 < jm) {
	      jm = j + jb - 1;
	    }
	    if (jm > jj) {
	      jm_m_jj = jm - jj;
	      dger_(&km,&jm_m_jj,&m_one,&ab[jj_ldab+kv+2],&ione,
		    &ab[kv+jj_ldab+ldab],&ldabm1,
		    &ab[kv+1+jj_ldab+ldab],&ldabm1);
	    } /* end if (jm > jj) */
	  } else {
		/*
               ELSE

                 If pivot is zero, set INFO to the index of the pivot
                 unless a zero pivot has already been found.
	         IF( INFO.EQ.0 )
     $               INFO = JJ
                 END IF
		*/
	    if (info == 0) {
	      info = jj;
	    }
	  } /* end else zero pivot */
	  /*	  
            Copy current column of A31 into the work array WORK31
               NW = MIN( JJ-J+1, I3 )
               IF( NW.GT.0 )
     $            CALL DCOPY( NW, AB( KV+KL+1-JJ+J, JJ ), 1,
     $                        WORK31( 1, JJ-J+1 ), 1 )
	  */
	  nw = i3;
	  if (jjmj +1 < nw) {
	    nw = jjmj + 1;
	  }
	  if (nw > 0) {
	    dcopy_(&nw,&ab[jj_ldab+kv+kl+1-jjmj],&ione,
		   &work31[1+(jjmj+1)*ldwork],&ione);
	  }
	  jj_ldab += ldab;
	} /* end for (jj...) */
	/*
   80       CONTINUE
            IF( J+JB.LE.N ) THEN
	*/
	if (j + jb <= n) {
	  /*
            Apply the row interchanges to the other blocks.
               J2 = MIN( JU-J+1, KV ) - JB
               J3 = MAX( 0, JU-J-KV+1 )
	  */
	  j2 = ju - j + 1;
	  if (kv < j2) {
	    j2 = kv;
	  }
	  j2 = j2 - jb;
	  j3 = 0;
	  if ((ju-j-kv+1) > j3) {
	    j3 = ju - j - kv+1;
	  }
	  /*
            Use DLASWP to apply the row interchanges to A12, A22, and
            A32.

               CALL DLASWP( J2, AB( KV+1-JB, J+JB ), LDAB-1, 1, JB,
     $                      IPIV( J ), 1 )
	  */
	  dlaswp_(&j2,&ab[kv+1-jb + jpjb_ldab],&ldabm1,&ione,&jb,
		  &ipiv[j],&ione);
	  /*
            Adjust the pivot indices.
*
               DO 90 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
   90          CONTINUE
	  */
	  for (i=j;i<=j+jb-1;i++) {
	    ipiv[i] = ipiv[i] + j - 1;
	  }
	  /*
            Apply the row interchanges to A13, A23, and A33
            columnwise.
               K2 = J - 1 + JB + J2
               DO 110 I = 1, J3
                  JJ = K2 + I
                  DO 100 II = J + I - 1, J + JB - 1
                     IP = IPIV( II )
                     IF( IP.NE.II ) THEN
                        TEMP = AB( KV+1+II-JJ, JJ )
                        AB( KV+1+II-JJ, JJ ) = AB( KV+1+IP-JJ, JJ )
                        AB( KV+1+IP-JJ, JJ ) = TEMP
                     END IF
  100             CONTINUE
  110          CONTINUE
	  */
	  k2 = j - 1 + jb + j2;
	  jj_ldab = k2 * ldab;
	  for (i=1;i<=j3;i++) {
	    jj = k2 + i;
	    for (ii = j + i -1 ;ii<=j+jb-1;ii++) {
	      ip = ipiv[ii];
	      if (ip != ii) {
		temp = ab[kv+1+ii-jj + jj_ldab];
		ab[kv+1+ii-jj+jj_ldab] = ab[kv+1+ip-jj + jj_ldab];
		ab[kv+1+ip-jj + jj_ldab] = temp;
	      }
	    } /* end for (ii..) */
	    jj_ldab += ldab;
	  } /* end for (i...) */
	  /*
            Update the relevant part of the trailing submatrix
               IF( J2.GT.0 ) THEN
	  */
	  if (j2 > 0) {
	    /*
              Update A12
	    CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit',
     $                        JB, J2, ONE, AB( KV+1, J ), LDAB-1,
     $                        AB( KV+1-JB, J+JB ), LDAB-1 )
	    */
	    dtrsm_(&l_char, &l_char, &n_char, &u_char,
		   &jb, &j2, &one, &ab[kv + 1 + j_ldab], &ldabm1,
		   &ab[kv+1-jb + jpjb_ldab], &ldabm1 );
/*
                  IF( I2.GT.0 ) THEN
*
*                    Update A22
*
	      CALL DGEMM( 'No transpose', 'No transpose', I2, J2,
     $                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1,
     $                           AB( KV+1-JB, J+JB ), LDAB-1, ONE,
     $                           AB( KV+1, J+JB ), LDAB-1 )
                  END IF
	      */
	    if (i2 > 0) {
	      /*
                Update A22
	      */
	      dgemm_(&n_char,&n_char,&i2,&j2,&jb,&m_one, 
		     &ab[j_ldab + kv +1 + jb],
		     &ldabm1,&ab[jpjb_ldab + kv+1 - jb], &ldabm1,
		     &one, &ab[jpjb_ldab + kv+1], &ldabm1);
	    }
	      /*
                  IF( I3.GT.0 ) THEN
*
*                    Update A32
*
                     CALL DGEMM( 'No transpose', 'No transpose', I3, J2,
     $                           JB, -ONE, WORK31, LDWORK,
     $                           AB( KV+1-JB, J+JB ), LDAB-1, ONE,
     $                           AB( KV+KL+1-JB, J+JB ), LDAB-1 )
                  END IF
	      */
	    if (i3 > 0){
	      dgemm_(&n_char, &n_char,&i3, &i2, &jb, &m_one, work31, &ldwork,
		     &ab[jpjb_ldab+kv+1-jb], &ldabm1, &one,
		     &ab[jpjb_ldab+kv+1-jb+kl], &ldabm1);

	    }
	    /*

               END IF
	    */
	  } /* end if j2 > 0 */
	  /*
               IF( J3.GT.0 ) THEN
	  */
	  if (j3 > 0) {
	    /*
	      Copy the lower triangle of A13 into the work array
              WORK13
                  DO 130 JJ = 1, J3
                     DO 120 II = JJ, JB
                        WORK13( II, JJ ) = AB( II-JJ+1, JJ+J+KV-1 )
  120                CONTINUE
  130             CONTINUE
	    */
	    jj_ldwork = ldwork;
	    jjpjpkv_ldab = kv*ldab + j_ldab;
	    for (jj=1;jj<=j3;jj++) {
	      for (ii= jj;ii<=jb;ii++) {
		work13[jj_ldwork + ii] = ab[jjpjpkv_ldab + ii-jj];
	      }
	      jj_ldwork += ldwork;
	      jjpjpkv_ldab += ldab;
	    }
	    /*
              Update A13 in the work array

                  CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit',
     $                        JB, J3, ONE, AB( KV+1, J ), LDAB-1,
     $                        WORK13, LDWORK )
	    */
	    dtrsm_(&l_char, &l_char, &n_char, &u_char,
		   &jb, &j3, &one, &ab[j_ldab + kv], 
		   &ldabm1,work13_p,&ldwork);
	    /*
*
                  IF( I2.GT.0 ) THEN
*
*                    Update A23
*
                     CALL DGEMM( 'No transpose', 'No transpose', I2, J3,
     $                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1,
     $                           WORK13, LDWORK, ONE, AB( 1+JB, J+KV ),
     $                           LDAB-1 )
                  END IF
	    */
	    jpkv_ldab = j_ldab + kv_ldab;
	    if (i2 > 0) {
	      /*
		Update A23
	      */
	      dgemm_(&n_char,&n_char,&i2,&j3,&jb,&m_one,
		     &ab[j_ldab + kv +jb], &ldabm1, work13_p,
		     &ldwork, &one, &ab[jpkv_ldab+jb],&ldabm1);
	    } /* end if (i2 > 0) */
	    /*

                  IF( I3.GT.0 ) THEN
*
*                    Update A33
*
                     CALL DGEMM( 'No transpose', 'No transpose', I3, J3,
     $                           JB, -ONE, WORK31, LDWORK, WORK13,
     $                           LDWORK, ONE, AB( 1+KL, J+KV ), LDAB-1 )
                  END IF
	    */
	    if (i3 > 0) {
	      /*
		Update A33
	      */
	      dgemm_ (&n_char,&n_char,&i3,&j3, &jb, &m_one, &work31_p[0],
		      &ldwork,work13_p,&ldwork,&one,
		      &ab[jpkv_ldab + kl + 1],&ldabm1);
	    }
	    /*
*                 Copy the lower triangle of A13 back into place
*
                  DO 150 JJ = 1, J3
                     DO 140 II = JJ, JB
                        AB( II-JJ+1, JJ+J+KV-1 ) = WORK13( II, JJ )
  140                CONTINUE
  150             CONTINUE
	    */
	    jj_ldab = ldab;
	    jj_ldwork = ldwork;
	    for (jj = 1;jj <- j3;jj ++ ) {
	      for (ii=jj;ii<=jb;ii++) {
		ab[jpkv_ldab + jj_ldab -ldab + ii - jj + 1] = work13[ii + jj_ldwork];
	      }
	      jj_ldab += ldab;
	      jj_ldwork += ldwork;
	    }
	    /*
               END IF
            ELSE
	    */
	  } /* end if (j3 > 0) */
	} else {
	  /*
	    j + jb >= n
	    Adjust the pivot indices.

               DO 160 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
  160          CONTINUE
            END IF
	  */
	  for (i=j;i<=j+jb-1;i++) {
	    ipiv[i] += (j - 1);
	  }
	} /* else (j + jb >= n) */
	/*
          Partially undo the interchanges in the current block to
          restore the upper triangular form of A31 and copy the upper
          triangle of A31 back into place
*
            DO 170 JJ = J + JB - 1, J, -1

               JP = IPIV( JJ ) - JJ + 1
               IF( JP.NE.1 ) THEN
	*/
	jj_ldab = jpjb_ldab - ldab;
	for (jj = j + jb - 1;jj >= j; jj--) {
	  jp = ipiv[jj] - jj + 1;
	  if (jp != 1) {
	    /*
              Apply interchange to columns J to JJ-1

                  IF( JP+JJ-1.LT.J+KL ) THEN
*
*                    The interchange does not affect A31
*
                     CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,
     $                           AB( KV+JP+JJ-J, J ), LDAB-1 )
                  ELSE
*
*                    The interchange does affect A31
*
                     CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,
     $                           WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                  END IF
               END IF
	    */
	    if ((jp + jj - 1) < (j + kl)) {
	      jjmj = jj - j;
	      dswap_(&jjmj,&ab[j_ldab+kv+1+jjmj],&ldabm1,
		     &ab[j_ldab+kv+jp+jjmj], &ldabm1);
	    } else {
	      /*
		The interchange does affect a31.
	      */
	      dswap_(&jjmj, &ab[j_ldab+kv+1+jjmj],&ldabm1,
		     &work31[ldwork+jp+jj-j-kl],&ldwork);
	    }
	  } /* end if (jp != 0) */
	  /*
            Copy the current column of A31 back into place

               NW = MIN( I3, JJ-J+1 )
               IF( NW.GT.0 )
     $            CALL DCOPY( NW, WORK31( 1, JJ-J+1 ), 1,
     $                        AB( KV+KL+1-JJ+J, JJ ), 1 )
	  */
	  nw = i3;
	  if ((jj-j + 1) < nw) {
	    nw = jj-j + 1;
	  }
	  if (nw > 0) {
	    dcopy_(&nw,&work31[1+(jj-j+1)*ldwork], &ione, 
		   &ab[jj_ldab + kv + kl + 1 + j - jj],&ione);
	  }
	  jj_ldab -= ldab;
	/*
  170       CONTINUE
	*/
	} /* end for jj */
	j_ldab += ldab;
	/*
  180    CONTINUE
	*/
      } /* end for(j...) */
    } /* end blocked mode */
  } /* end if valid input. */
  *infop = info;
}
