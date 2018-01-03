#include "system_includes.h"
#include "dgetrf.h"
#include "dgetrs.h"
#include "dtrsm.h"
#include "dlaswp.h"
#include "blas.h"

int main(int argc, char **argv)
{
  double result;
  double result2;
  double answer;
  double answer2;
  double two;
  double one;
  double zero;
  double neg_one;
  double msqrt;
  double dm;
  double *a;
  double *b;
  double *c;
  double cv[4];
  double *xsol;
  double *xtest;
  double *y;
  int64_t ask_for;
  int64_t one_l;
  int64_t col_m_o_2;
  int    *ipiv;
  int m;
  int msq;

  int mp1;
  int ldc;

  int iresult;
  int ianswer;

  int i;
  int j;

  int two_m_m1;
  int inc1;

  int itwo;
  int info;

  int mm1;
  int izero;

  char n_char;
  char l_char;
  char u_char;
  char r_char;
  
  m = 9;
  dm = (double)m;
  /*
    build m*m  matrix a
    off diagonal elements a[i,j] = -1 where i != j;
    a[i,i] = m+i, for i=0;m-1;
  */
  one_l   = (int64_t)1;
  zero    = 0.0;
  one     = 1.0;
  two     = 2.0;
  neg_one    = zero - one;
  inc1    = 1;
  itwo    = 2;
  izero   = 0;
  msq = m*m;
  mp1  = m + 1;
  two_m_m1 = m + m - 1;
  n_char = 'N';
  l_char = 'L';
  u_char = 'U';
  r_char = 'R';
  c = &cv[0];
  ask_for = sizeof(double) * msq;
  a  = (double*)calloc(one_l,ask_for);
  for (i=0;i<msq;i++) {
    a[i] = neg_one;
  }
  j = 0;
  for (i=0;i<msq;i+=mp1) {
    a[i] = m + j;
    j += 1;
  }
  /*
    build m*m  matrix b
    off diagonal elements b[i,j] = -1 where i != j;
    b[i,i] = 2m-1-i, for i=0;m-1;
  */
  ask_for = sizeof(double) * msq;
  b  = (double*)calloc(one_l,ask_for);
  for (i=0;i<msq;i++) {
    b[i] = neg_one;
  }
  j = 0;
  for (i=0;i<msq;i+=mp1) {
    b[i] = two_m_m1 - j;
    j += 1;
  }
  ask_for = sizeof(double) * m;
  xsol = (double*)calloc(one_l,ask_for);
  for (i=0;i<m;i++) {
    xsol[i] = one;
  }
  ask_for = sizeof(double) * m;
  xtest = (double*)calloc(one_l,ask_for);
  ask_for = sizeof(double) * m;
  y = (double*)calloc(one_l,ask_for);
  ask_for = sizeof(int) * m;
  ipiv = (int*)calloc(one_l,ask_for);
  /*
    test dnrm2.
  */
  msqrt = sqrt(dm);
  result   = dnrm2(&m, xsol, &inc1);
  fprintf(stdout," dnrm2 result should be %le, got %le\n",
	  msqrt,result);
  fflush(stdout);
  /*
    test dcopy.
  */
  dcopy(&m,xsol,&inc1,xtest,&inc1);
  /*
    Now xtest should be all 1's.
  */
  result = dnrm2(&m, xtest, &inc1);
  fprintf (stdout,"dcopy, nrm2 of xtest should be %le = %le\n",
	   msqrt,result);
  fflush(stdout);
  /*
    test dscal.
  */
  dscal(&m,&two,xtest,&inc1);
  /*
    Now xtest should be all 2's
  */
  answer = msqrt + msqrt;
  result = dnrm2(&m,xtest, &inc1);
  fprintf(stdout,"dscal test, nrm2 of xtest should be %le = %le\n",
	  answer,result);
  fflush(stdout);
  /*
    test daxpy. 
  */
  daxpy(&m,&two,xsol,&inc1,xtest,&inc1);
  /*
    now xtest should be all 4's
  */
  answer = answer + answer;
  result = dnrm2(&m,xtest, &inc1);
  fprintf(stdout,"daxpy test, nrm2 of xtest should be %le = %le\n",
	  answer,result);
  fflush(stdout);
  /*
    Test ddot, ddot(xtest,xsol) should be 4 * m
  */
  result = ddot(&m,xsol,&inc1,xtest,&inc1);
  answer = 4.0 * dm;
  fprintf(stdout,"ddot test,result should be %le = %le\n",
	  answer,result);
  fflush(stdout);
  /*
    Test idamax on col m/2 of a, result should be m/2
  */
  ianswer = m>>1;
  col_m_o_2 = ianswer * m;
  iresult = idamax(&m,&a[col_m_o_2],&inc1);
  fprintf(stdout,"idamax test on positive max abs, answer should be %d = %d\n",
	  ianswer,iresult);
  fflush(stdout);
  /*
    Now Test idamax on  - col m/2 of a, result should be m/2
  */
  dcopy(&m,&a[col_m_o_2],&inc1,xtest,&inc1);
  dscal(&m,&neg_one,xtest,&inc1);
  iresult = idamax(&m,xtest,&inc1);
  fprintf(stdout,"idamax test on negative max abs, answer should be %d = %d\n",
	  ianswer,iresult);
  fflush(stdout);
  /*
    Now test dgemv form y = a * xsol, should give a vector 
    with elements [1 .. m]
    first zero y.
  */
  for (i=0;i<m;i++) {
    y[i] = zero;
  }
  dgemv(&n_char,&m,&m,&one,a,&m,xsol,&inc1,&one,y,&inc1);
  fprintf(stdout,"dgemv test:\n answer\tresult\n");
  for (i=0;i<m;i++) {
    result = ((double)i) + one;
    fprintf(stdout,"%le\t%le\n",result,y[i]);
  }
  fflush(stdout);
  /*
    Now test dgemm on leading 2 by 2 submatrices of  a and b.
    Result should be
     2*m*m - m + 1 ,     -3m + 2
     -3m                  2*m*m - 1
  */
  ldc = 2;
  c[0] = 0;
  c[1] = 0;
  c[2] = 0;
  c[3] = 0;
  dgemm(&n_char,&n_char,&itwo,&itwo,&itwo,&one,a,&m,b,&m,&one,c,&ldc);
  answer = (two * (dm * dm)) - dm + one;
  answer2 = -3.0 * dm + two;
  
  fprintf(stdout, "dgemm test: Answer is \n"
	  "%le\t%le\n",answer,answer2);
  answer = answer2 - two;
  answer2 = (two * (dm *dm))  - one;
  fprintf(stdout, "%le\t%le\n",answer,answer2);
  fprintf(stdout, "\ndgemm test: result is \n"
	  "%le\t%le\n%le\t%le\n",c[0],c[2],c[1],c[3]);
  fflush(stdout);
  /*
    NB we only test DTRSM for the 'L' 'U' 'N' 'N',
                                  'L' 'U' 'N' 'U',
				  'L' 'L' 'N' 'N',
 				  'L' 'L' 'N' 'U',
    combinations. The other twelve really ought to be
                  tested as well.

  */
  /*
    Test dtrsm for nonunit upper  triangular solve.
    Upper triangular solve, non unit with a and solution of all ones.
    y will be odd numbers from 1 to 2m -1.
  */
  y[0] = one;
  for (i=1;i<m;i++) {
    y[i] = y[i-1] + two;
  }
  dcopy(&m,y,&inc1,xtest,&inc1);
  dtrsm(&l_char,&u_char,&n_char,&n_char,&m,&m,&one,a,&m,xtest,&m);
  daxpy(&m,&neg_one,xsol,&inc1,xtest,&inc1);
  answer = dnrm2(&m,xtest,&inc1);
  fprintf(stdout,"dtrsm test for nonunit uppper triangular solve\n");
  fprintf(stdout,"   nrm2(xtest - xsol ) = %le\n",answer);
  for (i=0;i<m;i++) {
    fprintf(stdout,"   xtest[%d] - xsol[%d] = %le\n",i,i,xtest[i]);
  }
  fflush(stdout);
  /*
    Test dtrsm for unit upper  triangular solve.
    Upper triangular solve, non unit with a and solution of all ones.
    y will be odd numbers from 1 to 2m -1.
  */
  y[0] = two - dm;
  for (i=1;i<m;i++) {
    y[i] = y[i-1] + one;
  }
  dcopy(&m,y,&inc1,xtest,&inc1);
  dtrsm(&l_char,&u_char,&n_char,&u_char,&m,&m,&one,a,&m,xtest,&m);
  daxpy(&m,&neg_one,xsol,&inc1,xtest,&inc1);
  answer = dnrm2(&m,xtest,&inc1);
  fprintf(stdout,"dtrsm test for unit uppper triangular solve\n");
  fprintf(stdout,"   nrm2(xtest - xsol ) = %le\n",answer);
  for (i=0;i<m;i++) {
    fprintf(stdout,"   xtest[%d] - xsol[%d] = %le\n",i,i,xtest[i]);
  }
  fflush(stdout);
  /*
    Test dtrsm for nonunit lower  triangular solve.
    Upper triangular solve, non unit with a and solution of all ones.
    y will be odd numbers from 1 to 2m -1.
  */
  for (i=0;i<m;i++) {
    y[i] = dm;
  }
  dcopy(&m,y,&inc1,xtest,&inc1);
  dtrsm(&l_char,&l_char,&n_char,&n_char,&m,&m,&one,a,&m,xtest,&m);
  daxpy(&m,&neg_one,xsol,&inc1,xtest,&inc1);
  answer = dnrm2(&m,xtest,&inc1);
  fprintf(stdout,"dtrsm test for nonunit lower triangular solve\n");
  fprintf(stdout,"   nrm2(xtest - xsol ) = %le\n",answer);
  for (i=0;i<m;i++) {
    fprintf(stdout,"   xtest[%d] - xsol[%d] = %le\n",i,i,xtest[i]);
  }
  fflush(stdout);
  /*
    Test dtrsm for unit lower  triangular solve.
    Upper triangular solve, non unit with a and solution of all ones.
    y will be odd numbers from 1 to 2m -1.
  */
  y[0] = one;
  for (i=1;i<m;i++) {
    y[i] = y[i-1] - one;
  }
  dcopy(&m,y,&inc1,xtest,&inc1);
  dtrsm(&l_char,&l_char,&n_char,&u_char,&m,&m,&one,a,&m,xtest,&m);
  daxpy(&m,&neg_one,xsol,&inc1,xtest,&inc1);
  answer = dnrm2(&m,xtest,&inc1);
  fprintf(stdout,"dtrsm test for nonunit lower triangular solve\n");
  fprintf(stdout,"   nrm2(xtest - xsol ) = %le\n",answer);
  for (i=0;i<m;i++) {
    fprintf(stdout,"   xtest[%d] - xsol[%d] = %le\n",i,i,xtest[i]);
  }
  fflush(stdout);
  /*
    Now test dgetrf and dgetrs with matrix a.
    Form y:
  */
  for (i=0;i<m;i++) {
    y[i] = zero;
  }
  dgemv(&n_char,&m,&m,&one,a,&m,xsol,&inc1,&one,y,&inc1);
  /*
    Factor. a.
  */
  dgetrf(&m,&m,a,&m,ipiv,&info);
  /*
    Check info code.
  */
  fprintf(stdout,"dgetrf on a, return code was %d\n",info);
  fflush(stdout);
  /*
    Now solve for xtest (should be = xsol).
  */
  dcopy(&m,y,&inc1,xtest,&inc1);
  dgetrs(&n_char,&m,&inc1,a,&m,ipiv,xtest,&m,&info);
  /*
    Check info code.
  */
  fprintf(stdout,"dgetrs on a, return code was %d\n",info);
  for (i=0;i<m;i++) {
    fprintf(stdout,"ipiv[%d] = %d\n",i,ipiv[i]);
  }
  fflush(stdout);
  daxpy(&m,&neg_one,xsol,&inc1,xtest,&inc1);
  answer = dnrm2(&m,xtest,&inc1);
  fprintf(stdout,"dgetrs test for a solve\n");
  fprintf(stdout,"   nrm2(xtest - xsol ) = %le\n",answer);
  for (i=0;i<m;i++) {
    fprintf(stdout,"   xtest[%d] - xsol[%d] = %le\n",i,i,xtest[i]);
  }
  fflush(stdout);
  
  /*
    Now test dgetrf and dgetrs with matrix b.
    Form y:
  */
  for (i=0;i<m;i++) {
    y[i] = zero;
  }
  dgemv(&n_char,&m,&m,&one,b,&m,xsol,&inc1,&one,y,&inc1);
  /*
    Factor. b.
  */
  dgetrf(&m,&m,b,&m,ipiv,&info);
  /*
    Check info code.
  */
  fprintf(stdout,"dgetrf on b, return code was %d\n",info);
  fflush(stdout);
  /*
    Now solve for xtest (should be = xsol).
  */
  dcopy(&m,y,&inc1,xtest,&inc1);
  dgetrs(&n_char,&m,&inc1,b,&m,ipiv,xtest,&m,&info);
  /*
    Check info code.
  */
  fprintf(stdout,"dgetrs on b, return code was %d\n",info);
  for (i=0;i<m;i++) {
    fprintf(stdout,"ipiv[%d] = %d\n",i,ipiv[i]);
  }
  fflush(stdout);
  daxpy(&m,&neg_one,xsol,&inc1,xtest,&inc1);
  answer = dnrm2(&m,xtest,&inc1);
  fprintf(stdout,"dgetrs test for b solve\n");
  fprintf(stdout,"   nrm2(xtest - xsol ) = %le\n",answer);
  for (i=0;i<m;i++) {
    fprintf(stdout,"   xtest[%d] - xsol[%d] = %le\n",i,i,xtest[i]);
  }
  fflush(stdout);
  /*
    Check to see if dlaswp works.
    First rebuild a.
  */
  for (i=0;i<msq;i++) {
    a[i] = neg_one;
  }
  j = 0;
  for (i=0;i<msq;i+=mp1) {
    a[i] = m + j;
    j += 1;
  }
  fprintf(stdout," a before dlaswp\n");
  fprintf(stdout,"%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t\n",
	  a[0],a[9],a[18],a[27],a[36],a[45],a[54],a[63],a[72]);
  fprintf(stdout,"%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t\n",
	  a[1],a[10],a[19],a[28],a[37],a[46],a[55],a[64],a[73]);
  fprintf(stdout,"%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t\n",
	  a[2],a[11],a[20],a[29],a[38],a[47],a[56],a[65],a[74]);
  fprintf(stdout,"%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t\n",
	  a[3],a[12],a[21],a[30],a[39],a[48],a[57],a[66],a[75]);
  fprintf(stdout,"%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t\n",
	  a[4],a[13],a[22],a[31],a[40],a[49],a[58],a[67],a[76]);
  fprintf(stdout,"%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t\n",
	  a[5],a[14],a[23],a[32],a[41],a[50],a[59],a[68],a[77]);
  fprintf(stdout,"%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t\n",
	  a[6],a[15],a[24],a[33],a[42],a[51],a[60],a[69],a[78]);
  fprintf(stdout,"%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t\n",
	  a[7],a[16],a[25],a[34],a[43],a[52],a[61],a[70],a[79]);
  fprintf(stdout,"%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t\n",
	  a[8],a[17],a[26],a[35],a[44],a[53],a[62],a[71],a[80]);
  ipiv[0] = 0;
  ipiv[1] = 1;
  ipiv[2] = 2;
  ipiv[3] = 5;
  ipiv[4] = 4;
  ipiv[5] = 5;
  ipiv[6] = 6;
  ipiv[7] = 7;
  mm1 = m -1;
  dlaswp(&m,a,&m,&izero,&mm1,ipiv,&inc1);
  fprintf(stdout," a after dlaswp\n");
  fprintf(stdout,"%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t\n",
	  a[0],a[9],a[18],a[27],a[36],a[45],a[54],a[63],a[72]);
  fprintf(stdout,"%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t\n",
	  a[1],a[10],a[19],a[28],a[37],a[46],a[55],a[64],a[73]);
  fprintf(stdout,"%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t\n",
	  a[2],a[11],a[20],a[29],a[38],a[47],a[56],a[65],a[74]);
  fprintf(stdout,"%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t\n",
	  a[3],a[12],a[21],a[30],a[39],a[48],a[57],a[66],a[75]);
  fprintf(stdout,"%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t\n",
	  a[4],a[13],a[22],a[31],a[40],a[49],a[58],a[67],a[76]);
  fprintf(stdout,"%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t\n",
	  a[5],a[14],a[23],a[32],a[41],a[50],a[59],a[68],a[77]);
  fprintf(stdout,"%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t\n",
	  a[6],a[15],a[24],a[33],a[42],a[51],a[60],a[69],a[78]);
  fprintf(stdout,"%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t\n",
	  a[7],a[16],a[25],a[34],a[43],a[52],a[61],a[70],a[79]);
  fprintf(stdout,"%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t\n",
	  a[8],a[17],a[26],a[35],a[44],a[53],a[62],a[71],a[80]);
  for (i=0;i<m;i++) {
    y[i] = zero;
  }
  dgemv(&n_char,&m,&m,&one,a,&m,xsol,&inc1,&one,y,&inc1);
  dcopy(&m,y,&inc1,xtest,&inc1);
  dgetrf(&m,&m,a,&m,ipiv,&info);
  /*
    Check info code.
  */
  fprintf(stdout,"dgetrf on p_a, return code was %d\n",info);
  fflush(stdout);
  /*
    Now solve for xtest (should be = xsol).
  */
  dgetrs(&n_char,&m,&inc1,a,&m,ipiv,xtest,&m,&info);
  /*
    Check info code.
  */
  fprintf(stdout,"dgetrs on p_a, return code was %d\n",info);
  for (i=0;i<m;i++) {
    fprintf(stdout,"ipiv[%d] = %d\n",i,ipiv[i]);
  }
  fflush(stdout);
  daxpy(&m,&neg_one,xsol,&inc1,xtest,&inc1);
  answer = dnrm2(&m,xtest,&inc1);
  fprintf(stdout,"dgetrs test forp_ a solve\n");
  fprintf(stdout,"   nrm2(xtest - xsol ) = %le\n",answer);
  for (i=0;i<m;i++) {
    fprintf(stdout,"   xtest[%d] - xsol[%d] = %le\n",i,i,xtest[i]);
  }
  fflush(stdout);
}
